function out = main_topo_tree(p)
% ============================================================
% ZigBee 多跳传输森林防火传感器网络仿真（函数化 + 路由可切换 + 指标输出）
% 输出 out 包含：
%   out.final: PDR/Delay/Hops/Alive/Energy/Lifetime 等最终指标
%   out.time : alive_ratio(t), delivered(t), generated(t) 等时间序列
%   out.diag : 诊断计数器
% ============================================================

% ---------- 默认参数 ----------
if nargin < 1, p = struct(); end
p = fill_defaults(p);

rng(p.seed);

%% ============ 生成拓扑（支持多种布点模式） ============
% 说明：本版本通过 generate_topology.m 支持 uniform/cluster/road 三种拓扑，
%      并可叠加 hotspot（“一个区域很多传感器”）
cfgTopo = struct();
cfgTopo.N = p.N;
cfgTopo.Lx = p.Lx;
cfgTopo.Ly = p.Ly;
cfgTopo.seed = p.seed;
cfgTopo.topo_mode = p.topo_mode;   % 'uniform' | 'cluster' | 'road'

% 基站/火区（用于热点中心默认值）
cfgTopo.BS_pos = p.BS_pos;
cfgTopo.fire_pos = p.fire_pos;
cfgTopo.Rf = p.Rf;

% hotspot（可选）
cfgTopo.hotspot_enable = p.hotspot_enable;
cfgTopo.hotspot_ratio  = p.hotspot_ratio;
cfgTopo.hotspot_sigma  = p.hotspot_sigma;
cfgTopo.hotspot_center = p.hotspot_center;

% cluster 参数
cfgTopo.Kc = p.Kc;
cfgTopo.cluster_sigma = p.cluster_sigma;
cfgTopo.cluster_center_mode = p.cluster_center_mode;
cfgTopo.cluster_mix_uniform_ratio = p.cluster_mix_uniform_ratio;

% road 参数（可选自定义 roads polyline）
cfgTopo.road_count = p.road_count;
cfgTopo.road_width = p.road_width;
cfgTopo.road_mix_uniform_ratio = p.road_mix_uniform_ratio;
cfgTopo.road_node_ratio = p.road_node_ratio;
if isfield(p,'roads') && ~isempty(p.roads)
    cfgTopo.roads = p.roads;
end

topo = generate_topology(cfgTopo);

pos = topo.pos_sensors;           % Nx2（仅传感器）
bs_idx = p.N + 1;
pos_all = topo.pos_all;           % (N+1)x2（最后一行为 BS）

% 邻接矩阵（距离阈值 Rc）
D = squareform(pdist(pos_all));
Adj = (D <= p.Rc) & (D > 0);

% 到基站距离（用于 greedy / mix）
d2bs = vecnorm(pos - p.BS_pos, 2, 2);
d2bs_all = [d2bs; 0];

% 火区节点标记
is_fire_node = vecnorm(pos - p.fire_pos, 2, 2) <= p.Rf;

%% ============ hop_to_bs 预计算（minhop 使用） ============
G = graph(Adj);
dist_hop = distances(G, 1:p.N, bs_idx);
hop_to_bs = inf(p.N+1,1);
hop_to_bs(1:p.N) = dist_hop(:);
hop_to_bs(bs_idx) = 0;

%% ============ 汇聚树（Converge Tree）预计算（tree 路由使用） ============
% 说明：
%   - 目的：把“汇聚网络”里常见的“面向基站的父子关系”显式建模出来。
%   - 做法：基于当前拓扑 + 期望链路代价（ETX 或 hop），从 BS 运行一次最短路，
%           得到每个节点到 BS 的父节点 parent(i)（即 i 的下一跳）。
%   - 容错：当 parent 失效/死亡/不相邻时，tree 路由会自动回退到 etx 选择。

% 期望 PRR（用于建树权重；与每次发送时的瞬时 prr 不同）
exp_prr = exp(-p.alpha * D);
exp_prr = max(exp_prr, p.prr_floor);
exp_prr = min(exp_prr, 1);

% 若用户仅设置 route_mode='tree_minhop'，则自动切换 tree_cost='hop'
if strcmpi(p.route_mode, 'tree_minhop')
    p.tree_cost = 'hop';
end

Adj_tree = Adj;
if p.USE_PRR_TH
    Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
end

parent = build_converge_tree(Adj_tree, exp_prr, bs_idx, p.tree_cost);
% parent(i)=0 表示到 BS 不连通（tree 模式会自动回退）

reachable_fire = isfinite(dist_hop) & is_fire_node;
if p.verbose
    fprintf("Sensors directly connected to BS: %d\n", sum(Adj(1:p.N, bs_idx)));
    fprintf("Connectivity check: fire sensors=%d, reachable=%d (%.1f%%)\n", ...
        sum(is_fire_node), sum(reachable_fire), 100*safe_div(sum(reachable_fire), sum(is_fire_node)));
    if any(isfinite(dist_hop))
        rh = dist_hop(isfinite(dist_hop));
        fprintf("Min hops (all->BS): %d, Median hops: %d, Max hops: %d\n", min(rh), median(rh), max(rh));
    end
end

%% ============ 状态初始化 ============
alive = true(p.N+1,1);
E = ones(p.N+1,1) * p.E0;
alive(bs_idx) = true; E(bs_idx) = inf;

Q = cell(p.N+1,1);
for i = 1:p.N+1, Q{i} = []; end

%% ============ 统计量 ============
gen_pkts   = 0;
deliv_pkts = 0;
drop_pkts  = 0;
sum_delay  = 0;
sum_hops   = 0;

diag.no_neighbor   = 0;
diag.ttl_drop      = 0;
diag.q_overflow    = 0;
diag.energy_dead   = 0;
diag.sent_cnt      = 0;
diag.succ_cnt      = 0;
diag.prr_drop      = 0;
diag.bs_in_nbr     = 0;
diag.try_to_bs     = 0;
diag.succ_to_bs    = 0;

% 用于能耗精确统计（不靠剩余能量反推）
energy.idle = 0;
energy.tx   = 0;
energy.rx   = 0;

% 时间序列（论文画寿命/随时间）
alive_ratio_ts = zeros(p.T,1);
gen_ts   = zeros(p.T,1);
deliv_ts = zeros(p.T,1);
drop_ts  = zeros(p.T,1);

% 寿命统计（FND/HND/LND）
FND = NaN; HND = NaN; LND = NaN;
dead_count_prev = 0;

%% ============ 主循环 ============
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

    % 3.5) 可选：周期性重建汇聚树（用于 tree 路由适配节点死亡）
    % 默认 def.tree_rebuild_period = inf（不重建，节省计算）
    if startsWith(lower(p.route_mode), 'tree') && isfinite(p.tree_rebuild_period)
        if mod(t-1, max(1,round(p.tree_rebuild_period))) == 0
            Adj_tree = Adj & (alive' & alive); % 只使用存活节点的链路
            if p.USE_PRR_TH
                Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
            end
            parent = build_converge_tree(Adj_tree, exp_prr, bs_idx, p.tree_cost);
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

    for s = 1:numel(senders)
        i = senders(s);
        if ~alive(i) || isempty(Q{i}), continue; end

        pkt = Q{i}(1); Q{i}(1) = [];

        % TTL
        pkt.ttl = pkt.ttl - 1;
        if pkt.ttl <= 0
            diag.ttl_drop = diag.ttl_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 邻居集合
        nbrs = find(Adj(i,:) & alive');
        if isempty(nbrs)
            diag.no_neighbor = diag.no_neighbor + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 计算 PRR
        d = D(i, nbrs);
        prr = exp(-p.alpha * d) + p.noise_sigma * randn(size(d));
        prr = max(prr, p.prr_floor);
        prr = min(prr, 1);

        % 链路阈值过滤
        if p.USE_PRR_TH
            keep = (prr >= p.PRR_MIN);
            nbrs = nbrs(keep);
            prr  = prr(keep);
            d    = d(keep);
            if isempty(nbrs)
                diag.no_neighbor = diag.no_neighbor + 1;
                drop_pkts = drop_pkts + 1;
                continue;
            end
        end

        diag.sent_cnt = diag.sent_cnt + 1;

        % 选择下一跳
        [j, prr_ij] = select_next_hop(p.route_mode, i, nbrs, prr, d, ...
                                     bs_idx, d2bs_all, E, hop_to_bs, parent, Adj, alive);

        if any(nbrs == bs_idx), diag.bs_in_nbr = diag.bs_in_nbr + 1; end
        if j == bs_idx, diag.try_to_bs = diag.try_to_bs + 1; end

        % 发送能耗（每次尝试都会扣 → 重传自动计入能耗）
        E(i) = E(i) - p.Etx;
        energy.tx = energy.tx + p.Etx;

        if E(i) <= 0 && alive(i)
            diag.energy_dead = diag.energy_dead + 1;
        end
        if E(i) <= 0
            alive(i) = false;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 成功/失败
        if rand < prr_ij
            diag.succ_cnt = diag.succ_cnt + 1;

            if j ~= bs_idx
                E(j) = E(j) - p.Erx;
                energy.rx = energy.rx + p.Erx;

                if E(j) <= 0 && alive(j)
                    diag.energy_dead = diag.energy_dead + 1;
                end
                if E(j) <= 0
                    alive(j) = false;
                end
            end

            pkt.hop = pkt.hop + 1;

            if j == bs_idx
                deliv_pkts = deliv_pkts + 1;
                diag.succ_to_bs = diag.succ_to_bs + 1;

                delay = t - pkt.born_t + 1;
                sum_delay = sum_delay + delay;
                sum_hops  = sum_hops  + pkt.hop;
            else
                [Q{j}, of] = push_pkt(Q{j}, pkt, p.Qmax);
                if of
                    diag.q_overflow = diag.q_overflow + 1;
                    drop_pkts = drop_pkts + 1;
                end
            end
        else
            diag.prr_drop = diag.prr_drop + 1;

            pkt.retx = pkt.retx + 1;
            if pkt.retx <= p.MAX_RETX
                % 更公平：建议放队尾（避免"卡死"占用发送机会）
                if p.RETX_TO_TAIL
                    Q{i} = [Q{i}, pkt];
                else
                    Q{i} = [pkt, Q{i}];
                end
            else
                drop_pkts = drop_pkts + 1;
            end
        end
    end

    % --------- 时间序列记录 ---------
    alive_ratio_ts(t) = mean(alive(1:p.N));
    gen_ts(t)   = gen_pkts;
    deliv_ts(t) = deliv_pkts;
    drop_ts(t)  = drop_pkts;

    % --------- 寿命统计 ---------
    dead_count = sum(~alive(1:p.N));
    if dead_count > dead_count_prev
        if isnan(FND) && dead_count >= 1
            FND = t;
        end
        if isnan(HND) && dead_count >= ceil(0.5*p.N)
            HND = t;
        end
        if isnan(LND) && dead_count >= p.N
            LND = t;
        end
        dead_count_prev = dead_count;
    end
end

%% ============ 汇总输出 ============
PDR = safe_div(deliv_pkts, gen_pkts);
avg_delay = safe_div(sum_delay, deliv_pkts);
avg_hops  = safe_div(sum_hops, deliv_pkts);
alive_ratio_final = mean(alive(1:p.N));

energy_total = energy.idle + energy.tx + energy.rx;

final = struct();
final.route_mode = p.route_mode;
final.PDR = PDR;
final.avg_delay = avg_delay;
final.avg_hops = avg_hops;
final.alive_ratio = alive_ratio_final;
final.generated = gen_pkts;
final.delivered = deliv_pkts;
final.dropped = drop_pkts;
final.energy_idle = energy.idle;
final.energy_tx = energy.tx;
final.energy_rx = energy.rx;
final.energy_total = energy_total;
final.FND = FND;
final.HND = HND;
final.LND = LND;

out = struct();
out.final = final;
out.time = struct("alive_ratio", alive_ratio_ts, "gen", gen_ts, "deliv", deliv_ts, "drop", drop_ts);
out.diag = diag;

if p.verbose
    fprintf("\n=== Simulation Results (route_mode=%s) ===\n", p.route_mode);
    fprintf("Generated packets: %d\n", gen_pkts);
    fprintf("Delivered packets:  %d\n", deliv_pkts);
    fprintf("Dropped packets:    %d\n", drop_pkts);
    fprintf("PDR: %.4f\n", PDR);
    fprintf("Avg delay (steps): %.2f\n", avg_delay);
    fprintf("Avg hops: %.2f\n", avg_hops);
    fprintf("Alive ratio (sensors): %.3f\n", alive_ratio_final);
    fprintf("Energy total (idle+tx+rx): %.2f (%.2f + %.2f + %.2f)\n", ...
        energy_total, energy.idle, energy.tx, energy.rx);
    fprintf("Lifetime: FND=%s, HND=%s, LND=%s\n", num2str(FND), num2str(HND), num2str(LND));
end

end

%% ==================== 路由选择（可切换策略） ====================
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
% parent(i)=0 表示不可达。
if startsWith(d_mode, "tree")
    j0 = parent(i);
    if j0 > 0 && alive(j0) && Adj(i,j0)
        j = j0;
        % 在候选 nbrs 里找对应 prr
        idx = find(nbrs == j, 1);
        if isempty(idx)
            % 理论上不会发生（因为 j0 要求 Adj(i,j0)=true 且 alive(j0)=true），做个最小保底
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
            [~, idx] = min(d2bs_n);
            j = nbrs(idx); prr_ij = prr(idx);
        end

    case 'etx'
        [~, idx] = min(etx);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'energy'
        score = 0.7 * normalize01(E_n) + 0.3 * normalize01(prr);
        [~, idx] = max(score);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'mindist'
        [~, idx] = min(d);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'minhop'
        [~, idx] = min(hop_n);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'mix'
        score = 0.7 * (-d2bs_n) + 0.3 * prr(:);
        [~, idx] = max(score);
        j = nbrs(idx); prr_ij = prr(idx);

    otherwise
        error("Unknown route_mode: %s", route_mode);
end
end

%% ==================== 入队（FIFO + 上限） ====================
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
    x = ones(size(v))*0.5;
else
    x = (v - mn) / (mx - mn);
end
end

%% ==================== 汇聚树构建（Dijkstra from BS） ====================
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
% 说明：
%   - 如果 parent(i) 不可用（不连通/死亡），会自动回退到 etx 选路，保证不崩。
def.tree_cost = 'etx';                 % 'etx' | 'hop'
def.tree_rebuild_period = inf;         % 每隔多少步重建一次（inf=不重建）
def.seed = 1;
def.verbose = false;

% ==================== 拓扑模式（2.1 补充） ====================
% 支持：'uniform'（随机均匀）/ 'cluster'（分簇）/ 'road'（沿道路）
def.topo_mode = 'uniform';

% hotspot：支持“一个区域很多传感器”（可叠加到任意 topo_mode）
def.hotspot_enable = false;
def.hotspot_ratio  = 0.25;                   % 热点节点比例
def.hotspot_sigma  = min(def.Lx,def.Ly)*0.05;% 热点扩散尺度（标准差）
def.hotspot_center = def.fire_pos;           % 默认以火区为热点中心

% cluster 参数
def.Kc = 6;                                  % 簇数量
def.cluster_sigma = min(def.Lx,def.Ly)*0.06; % 簇内扩散尺度
def.cluster_center_mode = 'random';          % 'random' 或 'grid'
def.cluster_mix_uniform_ratio = 0.2;         % 混入均匀节点比例

% road 参数
def.road_count = 2;                          % 自动生成道路数量（未提供 roads 时）
def.road_width = min(def.Lx,def.Ly)*0.01;    % 道路横向扩散尺度（标准差）
def.road_mix_uniform_ratio = 0.2;            % 混入均匀节点比例
def.road_node_ratio = 0.8;                   % 非均匀部分中，沿道路采样比例

% 重传放队尾（更推荐，避免"卡死"）
def.RETX_TO_TAIL = true;

% merge
fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(p, f)
        p.(f) = def.(f);
    end
end
end
