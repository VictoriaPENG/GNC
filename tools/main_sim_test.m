function main_sim_test()
% ============================================================
% ZigBee 多跳传输森林防火传感器网络：基础仿真（路由策略可切换版）
% 改动点：
% 1) 新增 route_mode 参数：可切换多种基线路由策略
% 2) 新增 select_next_hop()：集中实现下一跳选择逻辑
% 3) 新增 hop_to_bs 预计算：支持 minhop（最短跳数）策略
% ============================================================

clc; clear; close all;
rng(1);

%% ============ 参数区 ============
Lx = 1000; Ly = 1000;           % 区域大小 (m)
N  = 700;                       % 传感器数量
Rc = 150;                       % 通信半径 (m)

T = 2000;                       % 仿真步数
BS_pos   = [500, 500];          % 基站位置
fire_pos = [250, 250];          % 火源中心
Rf = 220;                       % 火区半径

lambda = 0.05;                  % 火区节点每步产生报警包的概率(0~1)

% 链路成功率模型：PRR = exp(-alpha*d) + noise
alpha = 0.01;                   % 路径衰减因子（越大 PRR 越差）
prr_floor = 0.08;               % PRR 下限（避免全 0）
noise_sigma = 0.02;             % PRR 噪声（森林遮挡/干扰）

% （可选）链路质量阈值：低于阈值的邻居不允许作为下一跳
USE_PRR_TH = false;
PRR_MIN = 0.10;                 % 例如 0.10；USE_PRR_TH=true 时生效

% 能耗模型
E0    = 100;                    % 初始能量
Etx   = 0.9;                    % 发送能耗
Erx   = 0.4;                    % 接收能耗
Eidle = 0.001;                  % 空闲能耗（每步）

% 节点失效模型
p_rand    = 1e-5;               % 随机失效概率/步
beta_fire = 0.05;               % 火烧毁指数因子（越大越容易烧毁）
fire_kill_scale = 1e-4;         % 火烧毁概率缩放（避免过猛）

% 队列与转发限制
Qmax   = 50;                    % 每个节点最大队列长度（溢出丢包）
TTLmax = 25;                    % 最大跳数（防环路/死循环）
Ksend = 20;                     % 每步最多选 Ksend 个发送者各转发 1 个包
MAX_RETX = 5;                   % 每一跳最多重传次数

% ================= 路由策略切换 =================
% 可选：
%   'greedy'  : 地理贪心（选择更接近基站的邻居）
%   'etx'     : ETX 最小（等价 PRR 最大）
%   'energy'  : 能量均衡（倾向选剩余能量高的邻居）
%   'mindist' : 最短距离（最近邻）
%   'minhop'  : 最短跳数（预计算 hop_to_bs，选 hop 最小邻居）
%   'mix'     : 你原来的混合打分（距离 + PRR）
route_mode = 'mix';
% =================================================


%% ============ 生成拓扑（基站作为 N+1 节点） ============
pos = [rand(N,1)*Lx, rand(N,1)*Ly];   % N 个传感器
bs_idx = N + 1;                       % 基站索引
pos_all = [pos; BS_pos];

% 距离矩阵与邻接矩阵（距离阈值）
D = squareform(pdist(pos_all));       % (N+1)x(N+1)
Adj = (D <= Rc) & (D > 0);

% 传感器到基站距离（路由打分可用）
d2bs = vecnorm(pos - BS_pos, 2, 2);
d2bs_all = [d2bs; 0];

% 火区节点（只对传感器 1..N）
is_fire_node = vecnorm(pos - fire_pos, 2, 2) <= Rf;

%% ============ 仿真开始前连通性检查（最短跳数） ============
% 用图论在"距离邻接图"上预计算每个节点到基站的 hop_to_bs（minhop 策略会用）
G = graph(Adj);
dist_hop = distances(G, 1:N, bs_idx);     % 传感器到基站的 hop 数（不可达为 Inf）
reachable = isfinite(dist_hop);
fprintf("Sensors directly connected to BS: %d\n", sum(Adj(1:N, bs_idx)));
fprintf("Connectivity check: fire sensors=%d, reachable=%d (%.1f%%)\n", ...
    sum(is_fire_node), sum(reachable & is_fire_node), ...
    100*safe_div(sum(reachable & is_fire_node), sum(is_fire_node)));

if any(reachable)
    rh = dist_hop(reachable);
    fprintf("Min hops (fire->BS): %d, Median hops: %d, Max hops: %d\n", ...
        min(rh), median(rh), max(rh));
end

hop_to_bs = inf(N+1,1);
hop_to_bs(1:N) = dist_hop(:);
hop_to_bs(bs_idx) = 0;

bs_neighbors = find(Adj(1:N, bs_idx));
if ~isempty(bs_neighbors)
    fprintf("Example BS-neighbor sensor index: %d\n", bs_neighbors(1));
end

%% ============ 状态初始化 ============
alive = true(N+1,1);
E = ones(N+1,1) * E0;

% 基站永远存活，能量无穷
alive(bs_idx) = true;
E(bs_idx) = inf;

% 队列：每个节点一个 cell，存储包结构体数组
Q = cell(N+1,1);
for i = 1:N+1
    Q{i} = [];
end

%% ============ 统计量 ============
gen_pkts   = 0;
deliv_pkts = 0;
drop_pkts  = 0;

sum_delay = 0;
sum_hops  = 0;

diag.no_neighbor   = 0;
diag.ttl_drop      = 0;
diag.q_overflow    = 0;
diag.energy_dead   = 0;
diag.sent_cnt      = 0;
diag.succ_cnt      = 0;
diag.prr_drop      = 0;
diag.prr_mean_sum  = 0;
diag.prr_min_sum   = 0;
diag.prr_max_sum   = 0;
diag.q_mean_sum    = 0;
diag.q_peak        = 0;

diag.bs_in_nbr  = 0;
diag.try_to_bs  = 0;
diag.succ_to_bs = 0;

%% ============ 主循环 ============
for t = 1:T

    % 1) 空闲能耗（只对传感器 1..N）
    E(1:N) = E(1:N) - Eidle;
    alive(1:N) = alive(1:N) & (E(1:N) > 0);

    % 2) 随机失效
    rand_fail = (rand(N,1) < p_rand);
    alive(1:N) = alive(1:N) & ~rand_fail;

    % 3) 火烧毁（只对火区且存活的传感器）
    fire_mask = is_fire_node & alive(1:N);
    if any(fire_mask)
        d_fire = vecnorm(pos(fire_mask,:) - fire_pos, 2, 2);
        p_fire = 1 - exp(-beta_fire .* (Rf - d_fire));
        p_fire = min(1, fire_kill_scale * p_fire);
        kill = rand(sum(fire_mask),1) < p_fire;
        idxs = find(fire_mask);
        alive(idxs(kill)) = false;
    end

    % 4) 火区节点产生报警包 -> 入队
    srcs = find(is_fire_node & alive(1:N) & (rand(N,1) < lambda));
    for k = 1:numel(srcs)
        i = srcs(k);
        gen_pkts = gen_pkts + 1;
        pkt = struct("born_t", t, "hop", 0, "ttl", TTLmax, "retx", 0, "src", i);
        [Q{i}, of] = push_pkt(Q{i}, pkt, Qmax);
        if of
            diag.q_overflow = diag.q_overflow + 1;
            drop_pkts = drop_pkts + 1;  % 产生即入队失败算丢
        end
    end

    % 5) 每步选择最多 Ksend 个"非空队列且存活"的节点发送 1 个包
    nonempty = find(alive & cellfun(@(x) ~isempty(x), Q));
    nonempty(nonempty == bs_idx) = []; % 基站不发
    if numel(nonempty) > Ksend
        senders = nonempty(randperm(numel(nonempty), Ksend));
    else
        senders = nonempty;
    end

    for s = 1:numel(senders)
        i = senders(s);

        if ~alive(i) || isempty(Q{i})
            continue;
        end

        % 取队首包（FIFO）
        pkt = Q{i}(1);
        Q{i}(1) = [];

        % TTL 检查
        pkt.ttl = pkt.ttl - 1;
        if pkt.ttl <= 0
            diag.ttl_drop = diag.ttl_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 邻居集合：距离邻接 + 存活
        nbrs = find(Adj(i,:) & alive');
        if isempty(nbrs)
            diag.no_neighbor = diag.no_neighbor + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 计算到邻居的 PRR（简化信道）
        d = D(i, nbrs);
        prr = exp(-alpha * d) + noise_sigma * randn(size(d));
        prr = max(prr, prr_floor);
        prr = min(prr, 1);

        % 可选：链路阈值过滤
        if USE_PRR_TH
            keep = (prr >= PRR_MIN);
            nbrs = nbrs(keep);
            prr  = prr(keep);
            d    = d(keep);
            if isempty(nbrs)
                diag.no_neighbor = diag.no_neighbor + 1;
                drop_pkts = drop_pkts + 1;
                continue;
            end
        end

        % PRR 统计
        diag.sent_cnt = diag.sent_cnt + 1;
        diag.prr_mean_sum = diag.prr_mean_sum + mean(prr);
        diag.prr_min_sum  = diag.prr_min_sum  + min(prr);
        diag.prr_max_sum  = diag.prr_max_sum  + max(prr);

        % 选择下一跳（可切换策略）
        [j, prr_ij] = select_next_hop(route_mode, i, nbrs, prr, d, ...
                                     bs_idx, d2bs_all, E, hop_to_bs);

        % 记录"基站在邻居中"的诊断
        if any(nbrs == bs_idx)
            diag.bs_in_nbr = diag.bs_in_nbr + 1;
        end
        if j == bs_idx
            diag.try_to_bs = diag.try_to_bs + 1;
        end

        % 发送能耗
        E(i) = E(i) - Etx;
        if E(i) <= 0 && alive(i)
            diag.energy_dead = diag.energy_dead + 1;
        end
        if E(i) <= 0
            alive(i) = false;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        % 链路成功/失败
        if rand < prr_ij
            % 接收成功
            diag.succ_cnt = diag.succ_cnt + 1;

            % 非基站扣接收能耗
            if j ~= bs_idx
                E(j) = E(j) - Erx;
                if E(j) <= 0 && alive(j)
                    diag.energy_dead = diag.energy_dead + 1;
                end
                if E(j) <= 0
                    alive(j) = false;
                end
            end

            % 跳数+1
            pkt.hop = pkt.hop + 1;

            if j == bs_idx
                deliv_pkts = deliv_pkts + 1;
                diag.succ_to_bs = diag.succ_to_bs + 1;

                delay = t - pkt.born_t + 1;
                sum_delay = sum_delay + delay;
                sum_hops  = sum_hops  + pkt.hop;
            else
                [Q{j}, of] = push_pkt(Q{j}, pkt, Qmax);
                if of
                    diag.q_overflow = diag.q_overflow + 1;
                    drop_pkts = drop_pkts + 1;
                end
            end
        else
            % 失败：MAC 重传
            diag.prr_drop = diag.prr_drop + 1;

            pkt.retx = pkt.retx + 1;
            if pkt.retx <= MAX_RETX
                % 放回队首（更像立即重传）
                Q{i} = [pkt, Q{i}];
            else
                drop_pkts = drop_pkts + 1;
            end
        end
    end

    % 每步统计队列情况
    q_len = cellfun(@numel, Q(1:N));
    diag.q_mean_sum = diag.q_mean_sum + mean(q_len);
    diag.q_peak = max(diag.q_peak, max(q_len));
end

%% ============ 输出结果 ============
PDR = safe_div(deliv_pkts, gen_pkts);
avg_delay = safe_div(sum_delay, deliv_pkts);
avg_hops  = safe_div(sum_hops, deliv_pkts);
alive_ratio = mean(alive(1:N));

fprintf("\n=== Simulation Results (route_mode=%s) ===\n", route_mode);
fprintf("Generated packets: %d\n", gen_pkts);
fprintf("Delivered packets:  %d\n", deliv_pkts);
fprintf("Dropped packets:    %d\n", drop_pkts);
fprintf("PDR: %.4f\n", PDR);
fprintf("Avg delay (steps): %.2f\n", avg_delay);
fprintf("Avg hops: %.2f\n", avg_hops);
fprintf("Alive ratio (sensors): %.3f\n", alive_ratio);

% ===== PRR 诊断输出 =====
if diag.sent_cnt > 0
    fprintf("\n=== Diagnostics (PRR) ===\n");
    fprintf("Tx attempts: %d, Rx success: %d, success ratio: %.3f\n", ...
        diag.sent_cnt, diag.succ_cnt, safe_div(diag.succ_cnt, diag.sent_cnt));
    fprintf("Mean PRR (avg over steps): %.3f\n", safe_div(diag.prr_mean_sum, diag.sent_cnt));
    fprintf("Min  PRR (avg over steps): %.3f\n", safe_div(diag.prr_min_sum,  diag.sent_cnt));
    fprintf("Max  PRR (avg over steps): %.3f\n", safe_div(diag.prr_max_sum,  diag.sent_cnt));
    fprintf("Link-fail events (incl retries): %d\n", diag.prr_drop);
end

fprintf("\n=== Diagnostics (Queue/TTL/Neighbor) ===\n");
fprintf("Queue overflow events: %d\n", diag.q_overflow);
fprintf("Avg queue length (time avg): %.2f\n", safe_div(diag.q_mean_sum, T));
fprintf("Peak queue length: %d (Qmax=%d)\n", diag.q_peak, Qmax);
fprintf("No-neighbor drops: %d\n", diag.no_neighbor);
fprintf("TTL drops: %d\n", diag.ttl_drop);
fprintf("Energy-dead events: %d\n", diag.energy_dead);

fprintf("\n=== Diagnostics (BS) ===\n");
fprintf("BS in neighbor count: %d\n", diag.bs_in_nbr);
fprintf("Try-to-BS count: %d\n", diag.try_to_bs);
fprintf("Succ-to-BS count: %d\n", diag.succ_to_bs);

end

%% ==================== 路由选择（可切换策略） ====================
function [j, prr_ij] = select_next_hop(route_mode, i, nbrs, prr, d, ...
                                      bs_idx, d2bs_all, E, hop_to_bs)
% 输入：
%   nbrs: 邻居索引列表
%   prr : 对应链路 PRR（向量）
%   d   : 对应距离 D(i,nbrs)（向量）
% 输出：
%   j: 选择的下一跳
%   prr_ij: 对应 PRR

% 规则 0：如果邻居里有基站，很多 ZigBee 汇聚网络会优先直达（可保底 PDR）
if any(nbrs == bs_idx)
    j = bs_idx;
    prr_ij = prr(nbrs == bs_idx);
    return;
end

% 常用派生量
etx = 1 ./ max(prr, 1e-6);           % ETX ~ 1/PRR
d2bs_n = d2bs_all(nbrs);             % 邻居到基站距离
E_n = E(nbrs);                       % 邻居剩余能量
hop_n = hop_to_bs(nbrs);             % 邻居到基站 hop（预计算）

switch lower(route_mode)

    case 'greedy'
        % 地理贪心：选更接近基站的邻居（避免走远）
        % 若无更近邻居，就退化为选择 d2bs 最小（仍可推进）
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
        % ETX 最小（PRR 最大）
        [~, idx] = min(etx);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'energy'
        % 能量均衡：倾向选能量更高的邻居
        % 可在这里叠加可靠性：例如 score = 0.7*E + 0.3*PRR
        score = 0.7 * normalize01(E_n) + 0.3 * normalize01(prr);
        [~, idx] = max(score);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'mindist'
        % 最短距离（最近邻）
        [~, idx] = min(d);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'minhop'
        % 最短跳数：选 hop_to_bs 最小的邻居
        % 若存在不可达(Inf)，自然会排到后面
        [~, idx] = min(hop_n);
        j = nbrs(idx); prr_ij = prr(idx);

    case 'mix'
        % 你原来的混合打分：越靠近基站越好 + PRR 越高越好
        score = 0.7 * (-d2bs_n) + 0.3 * prr(:);
        [~, idx] = max(score);
        j = nbrs(idx); prr_ij = prr(idx);

    otherwise
        error("Unknown route_mode: %s", route_mode);
end
end

%% ==================== 队列入队（FIFO + 上限） ====================
function [Qnew, overflow] = push_pkt(Qold, pkt, Qmax)
overflow = false;
Qnew = Qold;
if numel(Qnew) >= Qmax
    overflow = true;
    return;
end
Qnew = [Qnew, pkt];
end

%% ==================== 安全除法 ====================
function y = safe_div(a,b)
if b == 0
    y = 0;
else
    y = a / b;
end
end

%% ==================== 归一化到 [0,1]（用于加权打分） ====================
function x = normalize01(v)
v = v(:);
mn = min(v); mx = max(v);
if mx - mn < 1e-12
    x = ones(size(v))*0.5;
else
    x = (v - mn) / (mx - mn);
end
end
