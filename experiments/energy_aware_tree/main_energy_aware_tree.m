function out = main_energy_aware_tree(p)
% ============================================================
% ZigBee 多跳传输森林防火传感器网络仿真（函数化 + 路由可切换 + 指标输出）
% 输出 out 包含：
%   out.final: PDR/Delay/Hops/Alive/Energy/Lifetime 等最终指标
%   out.time : alive_ratio(t), delivered(t), generated(t) 等时间序列
%   out.diag : 诊断计数器
%
% [新增] 能量感知建树（Energy-aware Tree）：
%   - 当 route_mode 为 'tree_energy' 或 p.tree_energy_aware=true 时启用
%   - 建树边权：ETX/hop * 能量惩罚(父节点残能越低代价越大)
%   - 建议配合 tree_rebuild_period 周期性重建
% ============================================================

% ---------- 默认参数 ----------
if nargin < 1, p = struct(); end
p = fill_defaults(p);

rng(p.seed);


%% ================= 拓扑生成（适配你现有 generate_topology） =================
topo = generate_topology(p);

% 所有节点坐标（N 个传感器 + 1 个 BS，BS 在最后）
pos = topo.pos_all;
bs_idx = size(pos, 1);

% 传感器坐标
sensor_pos = topo.pos_sensors;

% 火区节点判定（只针对传感器 1..N）
is_fire_node = (vecnorm(sensor_pos - p.fire_pos, 2, 2) <= p.Rf);

% 若后面代码需要与 pos 对齐，可用：
% is_fire_node_all = [is_fire_node; false];


% 节点编号：1..N 为传感器，bs_idx=N+1 为基站
p.N = size(pos,1)-1;
bs_idx = p.N + 1;

% 位置
BS_pos = pos(bs_idx,:);

%% ============ 邻接/距离 ============
D = squareform(pdist(pos));
Adj = (D > 0) & (D <= p.Rc);

% 基站不需要发包队列，但依然在图里参与连通性/路由
Adj(bs_idx, bs_idx) = false;

% ============ 预计算到 BS 的几何距离 ============
d2bs_all = vecnorm(pos - BS_pos, 2, 2);

%% ============ 预计算 hop 到 BS（用于诊断/部分路由） ============
hop_to_bs = inf(p.N+1,1);
hop_to_bs(bs_idx) = 0;
% BFS
q = bs_idx;
visited = false(p.N+1,1); visited(bs_idx)=true;
while ~isempty(q)
    u = q(1); q(1)=[];
    nbrs = find(Adj(u,:));
    for v = nbrs
        if ~visited(v)
            visited(v)=true;
            hop_to_bs(v)=hop_to_bs(u)+1;
            q(end+1)=v; %#ok<AGROW>
        end
    end
end

%% ============ 汇聚树（Converge Tree）预计算（tree 路由使用） ============
% 期望 PRR（用于建树权重；与每次发送时的瞬时 prr 不同）
exp_prr = exp(-p.alpha * D);
exp_prr = max(exp_prr, p.prr_floor);
exp_prr = min(exp_prr, 1);

% 若设置 minhop 类 tree，自动切换 tree_cost='hop'
if contains(lower(string(p.route_mode)), "minhop")
    p.tree_cost = 'hop';
end

Adj_tree = Adj;
if p.USE_PRR_TH
    Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
end

% 注意：此时 E 尚未进入仿真循环，先用 E0 初始化一次
E_init = p.E0 * ones(p.N+1,1);

if p.tree_energy_aware || contains(lower(string(p.route_mode)), "energy")
    parent = build_converge_tree_energy(Adj_tree, exp_prr, bs_idx, p.tree_cost, ...
        E_init, p.E0, p.tree_energy_beta, p.tree_energy_min_frac);
else
    parent = build_converge_tree(Adj_tree, exp_prr, bs_idx, p.tree_cost);
end
% parent(i)=0 表示到 BS 不连通（tree 模式会自动回退）

%% ============ 状态初始化 ============
alive = true(p.N+1,1);
E = ones(p.N+1,1) * p.E0;

% 队列（每个节点一个 FIFO）
Q = cell(p.N+1,1);
for i = 1:p.N+1, Q{i} = []; end

% 统计量
gen_pkts = 0;
del_pkts = 0;
drop_pkts = 0;

sum_delay = 0;
sum_hops = 0;

energy = struct('tx',0,'rx',0,'idle',0);
diag = struct('q_overflow',0,'ttl_drop',0,'dead_drop',0,'no_route_drop',0,'prr_fail',0,'retx_exceed',0);

% 时间序列输出
time.alive_ratio = zeros(p.T,1);
time.generated = zeros(p.T,1);
time.delivered = zeros(p.T,1);

%% ============ 仿真主循环 ============
for t = 1:p.T

    % 1) 空闲能耗
    E(1:p.N) = E(1:p.N) - p.Eidle;
    energy.idle = energy.idle + p.N * p.Eidle;

    alive(1:p.N) = alive(1:p.N) & (E(1:p.N) > 0);

    % 2) 随机失效
    rand_fail = (rand(p.N,1) < p.p_rand);
    alive(1:p.N) = alive(1:p.N) & ~rand_fail;

    % 3) 火烧毁
    fire_mask = is_fire_node & alive(1:p.N);
    if any(fire_mask)
        d_fire = vecnorm(pos(fire_mask,:) - p.fire_pos, 2, 2);
        p_fire = 1 - exp(-p.beta_fire .* (p.Rf - d_fire));
        p_fire = min(1, p.fire_kill_scale * p_fire);
        kill = rand(sum(fire_mask),1) < p_fire;
        idxs = find(fire_mask);
        alive(idxs(kill)) = false;
    end

    % 3.5) 周期性重建汇聚树（适配节点死亡/能量变化）
    if startsWith(lower(p.route_mode), 'tree') && isfinite(p.tree_rebuild_period)
        if mod(t-1, max(1,round(p.tree_rebuild_period))) == 0
            Adj_tree = Adj & (alive' & alive); % 只使用存活节点链路
            if p.USE_PRR_TH
                Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
            end

            if p.tree_energy_aware || contains(lower(string(p.route_mode)), "energy")
                parent = build_converge_tree_energy(Adj_tree, exp_prr, bs_idx, p.tree_cost, ...
                    E, p.E0, p.tree_energy_beta, p.tree_energy_min_frac);
            else
                parent = build_converge_tree(Adj_tree, exp_prr, bs_idx, p.tree_cost);
            end
        end
    end

    % 4) 产生报警包 -> 入队
    srcs = find(is_fire_node & alive(1:p.N) & (rand(p.N,1) < p.lambda));
    for k = 1:numel(srcs)
        i = srcs(k);
        gen_pkts = gen_pkts + 1;
        pkt = struct("born_t", t, "hop", 0, "ttl", p.TTLmax, "retx", 0, "src", i);
        [Q{i}, of] = push_pkt(Q{i}, pkt, p.Qmax);
        if of
            diag.q_overflow = diag.q_overflow + 1;
            drop_pkts = drop_pkts + 1;
        end
    end

    % 5) 每步选择最多 Ksend 个发送者
    nonempty = find(alive & cellfun(@(x) ~isempty(x), Q));
    nonempty(nonempty == bs_idx) = [];
    if numel(nonempty) > p.Ksend
        senders = nonempty(randperm(numel(nonempty), p.Ksend));
    else
        senders = nonempty;
    end

    % 6) 逐个发送（简化 MAC：每步每节点最多发 1 个包）
    for s = 1:numel(senders)
        i = senders(s);

        if isempty(Q{i}) || ~alive(i)
            continue;
        end

        pkt = Q{i}(1);
        Q{i}(1) = [];

        % TTL
        pkt.ttl = pkt.ttl - 1;
        if pkt.ttl <= 0
            diag.ttl_drop = diag.ttl_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 邻居集合（仅存活）
        nbrs = find(Adj(i,:) & alive');
        if isempty(nbrs)
            diag.no_route_drop = diag.no_route_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 瞬时 PRR：期望 PRR + 噪声
        prr = exp_prr(i, nbrs) + p.noise_sigma * randn(size(nbrs));
        prr = min(1, max(p.prr_floor, prr));

        % 选下一跳
        [j, prr_ij] = select_next_hop(p.route_mode, i, nbrs, prr, D(i,nbrs), ...
            bs_idx, d2bs_all, E, hop_to_bs, parent, Adj, alive);

        if j == 0
            diag.no_route_drop = diag.no_route_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 发送一次（含重传）
        success = false;
        for rtx = 0:p.MAX_RETX
            % TX/RX 能耗
            E(i) = E(i) - p.Etx;
            energy.tx = energy.tx + p.Etx;
            if j ~= bs_idx
                E(j) = E(j) - p.Erx;
                energy.rx = energy.rx + p.Erx;
            end

            % 发送后可能死亡
            if E(i) <= 0, alive(i) = false; end
            if j ~= bs_idx && E(j) <= 0, alive(j) = false; end

            if ~alive(i)
                diag.dead_drop = diag.dead_drop + 1;
                drop_pkts = drop_pkts + 1;
                break;
            end
            if j ~= bs_idx && ~alive(j)
                diag.dead_drop = diag.dead_drop + 1;
                drop_pkts = drop_pkts + 1;
                break;
            end

            % 链路成功/失败
            if rand() < prr_ij
                success = true;
                break;
            else
                diag.prr_fail = diag.prr_fail + 1;
                pkt.retx = pkt.retx + 1;
                if pkt.retx > p.MAX_RETX
                    diag.retx_exceed = diag.retx_exceed + 1;
                    drop_pkts = drop_pkts + 1;
                    break;
                end
            end
        end

        if ~success
            continue;
        end

        % 成功：更新 hop
        pkt.hop = pkt.hop + 1;

        % 到达 BS
        if j == bs_idx
            del_pkts = del_pkts + 1;
            sum_delay = sum_delay + (t - pkt.born_t);
            sum_hops = sum_hops + pkt.hop;
        else
            % 中继入队
            [Q{j}, of] = push_pkt(Q{j}, pkt, p.Qmax);
            if of
                diag.q_overflow = diag.q_overflow + 1;
                drop_pkts = drop_pkts + 1;
            end
        end
    end

    % 记录时间序列
    time.alive_ratio(t) = mean(alive(1:p.N));
    time.generated(t) = gen_pkts;
    time.delivered(t) = del_pkts;
end

%% ============ 输出汇总 ============
final = struct();
final.PDR = safe_div(del_pkts, gen_pkts);
final.Delay = safe_div(sum_delay, del_pkts);
final.Hops = safe_div(sum_hops, del_pkts);
final.Alive = mean(alive(1:p.N));
final.Eavg = mean(max(E(1:p.N),0));
final.Esum = sum(max(E(1:p.N),0));

% Lifetime（第一次死亡 / 半数死亡 / 全部死亡）简单估计
alive_ratio = time.alive_ratio;
fnd = find(alive_ratio < (1 - 1/p.N), 1, 'first');
hnd = find(alive_ratio <= 0.5, 1, 'first');
lnd = find(alive_ratio <= 0.0, 1, 'first');
if isempty(fnd), fnd = p.T; end
if isempty(hnd), hnd = p.T; end
if isempty(lnd), lnd = p.T; end
final.FND = fnd;
final.HND = hnd;
final.LND = lnd;

out = struct();
out.final = final;
out.time = time;
out.diag = diag;
out.energy = energy;
out.parent = parent;
out.param = p;

end

%% ==================== 选路函数 ====================
function [j, prr_ij] = select_next_hop(route_mode, i, nbrs, prr, d, ...
                                      bs_idx, d2bs_all, E, hop_to_bs, parent, Adj, alive)
% 基站优先保底（符合汇聚网络常见实现）
if any(nbrs == bs_idx)
    j = bs_idx;
    prr_ij = prr(nbrs == bs_idx);
    return;
end

d_mode = lower(string(route_mode));
etx = 1 ./ max(prr, 1e-6);
d2bs_n = d2bs_all(nbrs);
E_n = E(nbrs);
hop_n = hop_to_bs(nbrs);

% ---------------- tree（汇聚树）模式 ----------------
% route_mode:
%   - 'tree' / 'tree_etx' : 走 parent(i)（树由 ETX 建）
%   - 'tree_minhop'       : 走 parent(i)（树由 hop 建）
%   - 'tree_energy'       : 走 parent(i)（树由 ETX + 能量惩罚 建）
% parent(i)=0 表示不可达。
if startsWith(d_mode, "tree")
    j0 = parent(i);
    if j0 > 0 && alive(j0) && Adj(i,j0)
        j = j0;
        idx = find(nbrs == j, 1);
        if isempty(idx)
            prr_ij = min(prr);
        else
            prr_ij = prr(idx);
        end
        return;
    end
    % parent 不可用则回退到 etx（更鲁棒）
    d_mode = "etx";
end

switch char(d_mode)
    case 'greedy'
        cur = d2bs_all(i);
        better = (d2bs_n < cur);
        if any(better)
            [~, idx] = min(d2bs_n(better));
            cand = nbrs(better);
            j = cand(idx);
            prr_ij = prr(better); prr_ij = prr_ij(idx);
        else
            j = 0; prr_ij = 0;
        end

    case 'minhop'
        [~, idx] = min(hop_n);
        j = nbrs(idx);
        prr_ij = prr(idx);

    case 'mindist'
        [~, idx] = min(d);
        j = nbrs(idx);
        prr_ij = prr(idx);

    case 'energy'
        % 传统能量感知：偏向高能量 + 好链路
        score = normalize01(E_n) + normalize01(prr);
        [~, idx] = max(score);
        j = nbrs(idx);
        prr_ij = prr(idx);

    case 'mix'
        % 混合：好链路 + 离 BS 近 + 高能量
        score = 0.45*normalize01(prr) + 0.35*(1-normalize01(d2bs_n)) + 0.20*normalize01(E_n);
        [~, idx] = max(score);
        j = nbrs(idx);
        prr_ij = prr(idx);

    otherwise % 'etx'
        [~, idx] = min(etx);
        j = nbrs(idx);
        prr_ij = prr(idx);
end
end

%% ==================== 入队 ====================
function [Qnew, overflow] = push_pkt(Qold, pkt, Qmax)
overflow = false;
Qnew = Qold;
if numel(Qnew) >= Qmax
    overflow = true;
    return;
end
Qnew = [Qnew, pkt];
end

function y = safe_div(a,b)
if b == 0, y = 0; else, y = a / b; end
end

function x = normalize01(v)
v = v(:);
mn = min(v); mx = max(v);
if mx - mn < 1e-12
    x = zeros(size(v));
else
    x = (v - mn) / (mx - mn);
end
end

%% ==================== 汇聚树构建（原版） ====================
function parent = build_converge_tree(Adj, exp_prr, bs_idx, tree_cost)
% parent(i) 为节点 i 的下一跳（指向 BS 方向的父节点）。
% 不连通则 parent(i)=0。

n = size(Adj,1);
parent = zeros(n,1);

tree_cost = lower(string(tree_cost));
use_hop = (tree_cost == "hop") || (tree_cost == "minhop");

% Dijkstra (O(n^2)；n~700 足够快，且避免版本差异导致的 graph API 不兼容)
dist = inf(n,1);
visited = false(n,1);
dist(bs_idx) = 0;

for iter = 1:n
    % pick unvisited min dist
    du = dist;
    du(visited) = inf;
    [~, u] = min(du);
    if ~isfinite(dist(u)), break; end
    visited(u) = true;

    nbrs = find(Adj(u,:));
    if isempty(nbrs), continue; end

    if use_hop
        w = ones(size(nbrs));
    else
        p_ = exp_prr(u, nbrs);
        w = 1 ./ max(p_, 1e-6); % ETX
    end

    alt = dist(u) + w(:);
    better = alt < dist(nbrs);
    if any(better)
        idxs = nbrs(better);
        dist(idxs) = alt(better);
        parent(idxs) = u;
    end
end

% bs 自己不需要父节点
parent(bs_idx) = 0;
end

%% ==================== [新增] 能量感知汇聚树构建 ====================
function parent = build_converge_tree_energy(Adj, exp_prr, bs_idx, tree_cost, E, E0, beta, min_frac)
% 能量感知汇聚树（Energy-aware Tree）
% 思路：仍以"到 BS 的最短路"为主，但在边权中加入"父节点能量惩罚"，
%       使得更倾向选择残余能量更高的父节点作为下一跳（从而均衡转发负载）。
%
% 代价设计（u 作为父节点时）：
%   w(u->v) = base_w(u->v) * (1 + beta * (1 - clamp(E(u)/E0, min_frac, 1)))
%   - beta 越大越强调能量均衡
%   - min_frac 防止 E 很小导致权重过大而数值不稳
%
% 说明：这里只使用"当前残余能量 E"做静态权重；可配合 tree_rebuild_period 周期性重建以适配能量变化。

n = size(Adj,1);
parent = zeros(n,1);

tree_cost = lower(string(tree_cost));
use_hop = (tree_cost == "hop") || (tree_cost == "minhop");

% Dijkstra (O(n^2))
dist = inf(n,1);
visited = false(n,1);
dist(bs_idx) = 0;

% 安全处理
if nargin < 7 || isempty(beta), beta = 2.0; end
if nargin < 8 || isempty(min_frac), min_frac = 0.10; end
if isempty(E), E = E0 * ones(n,1); end

for iter = 1:n
    du = dist;
    du(visited) = inf;
    [~, u] = min(du);
    if ~isfinite(dist(u)), break; end
    visited(u) = true;

    nbrs = find(Adj(u,:));
    if isempty(nbrs), continue; end

    if use_hop
        base_w = ones(size(nbrs));
    else
        p_ = exp_prr(u, nbrs);
        base_w = 1 ./ max(p_, 1e-6); % ETX
    end

    % 父节点 u 的能量惩罚（BS 不惩罚）
    if u == bs_idx
        pen = 1.0;
    else
        frac = max(0, E(u)) / max(E0, 1e-9);
        frac = min(1.0, max(min_frac, frac));
        pen = 1.0 + beta * (1.0 - frac);
    end

    w = base_w * pen;

    alt = dist(u) + w(:);
    better = alt < dist(nbrs);
    if any(better)
        idxs = nbrs(better);
        dist(idxs) = alt(better);
        parent(idxs) = u;
    end
end

parent(bs_idx) = 0;
end

%% ==================== 默认参数填充 ====================
function p = fill_defaults(p)
def = struct();
def.Lx = 1000; def.Ly = 1000;
def.N  = 700;
def.Rc = 150;
def.T  = 2000;

def.BS_pos   = [500, 500];
def.fire_pos = [250, 250];
def.Rf = 220;

def.lambda = 0.05;

def.alpha = 0.01;
def.prr_floor = 0.08;
def.noise_sigma = 0.02;

def.USE_PRR_TH = true;
def.PRR_MIN = 0.15;

def.E0    = 100;
def.Etx   = 0.9;
def.Erx   = 0.4;
def.Eidle = 0.001;

def.p_rand = 1e-5;
def.beta_fire = 0.05;
def.fire_kill_scale = 1e-4;

def.Qmax   = 50;
def.TTLmax = 25;
def.Ksend = 20;
def.MAX_RETX = 5;

def.route_mode = 'mix';

% ==================== 汇聚树（Converge Tree）路由参数 ====================
% route_mode 取值：
%   - 'tree' / 'tree_etx' : 使用 ETX 代价构建汇聚树，转发固定走 parent(i)
%   - 'tree_minhop'       : 使用 hop 代价构建汇聚树，转发固定走 parent(i)
%   - 'tree_energy'       : 使用 ETX + 能量惩罚构建汇聚树（Energy-aware Tree）
% 说明：
%   - 如果 parent(i) 不可用（不连通/死亡），会自动回退到 etx 选路，保证不崩。
def.tree_cost = 'etx';                 % 'etx' | 'hop'
def.tree_rebuild_period = inf;         % 每隔多少步重建一次（inf=不重建）
def.tree_energy_aware = false;         % 是否启用能量感知建树（Energy-aware Tree）
def.tree_energy_beta  = 2.0;           % 能量惩罚强度（越大越偏向选择高能量父节点）
def.tree_energy_min_frac = 0.10;       % 残余能量归一化下限（避免权重过大）
def.seed = 1;
def.verbose = false;

% 合并默认
fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(p, f) || isempty(p.(f))
        p.(f) = def.(f);
    end
end
end
