function out = main_aodv_energy_hybrid(p)
% ============================================================
% AODV-local-switch + Dynamic Energy-aware Tree Hybrid Routing
% ------------------------------------------------------------
% Simplifies AODV RREQ/RREP into a local parent switching rule:
%   - Each node keeps a parent toward BS (from energy-aware tree)
%   - If parent becomes weak (PRR/energy) or unavailable, switch
%     to a better neighbor that is closer to BS (hop/distance).
%   - Periodic tree rebuild keeps global structure dynamic.
% ============================================================

% ---------- Default parameters ----------
if nargin < 1, p = struct(); end
p = fill_defaults(p);

% Route-mode based auto-enable
if contains(lower(string(p.route_mode)), "aodv")
    p.local_switch_enable = true;
end

rng(p.seed);

%% ================= Topology =================
topo = generate_topology(p);

pos = topo.pos_all;
bs_idx = size(pos, 1);
sensor_pos = topo.pos_sensors;

is_fire_node = (vecnorm(sensor_pos - p.fire_pos, 2, 2) <= p.Rf);

p.N = size(pos,1) - 1;
bs_idx = p.N + 1;
BS_pos = pos(bs_idx,:);

%% ============ Adjacency/Distance ============
D = squareform(pdist(pos));
Adj = (D > 0) & (D <= p.Rc);
Adj(bs_idx, bs_idx) = false;

%% ============ Precompute to BS ============
d2bs_all = vecnorm(pos - BS_pos, 2, 2);

%% ============ Precompute hop to BS ============
hop_to_bs = inf(p.N+1,1);
hop_to_bs(bs_idx) = 0;
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

%% ============ Converge Tree (Energy-aware) ============
exp_prr = exp(-p.alpha * D);
exp_prr = max(exp_prr, p.prr_floor);
exp_prr = min(exp_prr, 1);

if contains(lower(string(p.route_mode)), "minhop")
    p.tree_cost = 'hop';
end

Adj_tree = Adj;
if p.USE_PRR_TH
    Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
end

E_init = p.E0 * ones(p.N+1,1);
parent = build_converge_tree_energy(Adj_tree, exp_prr, bs_idx, p.tree_cost, ...
    E_init, p.E0, p.tree_energy_beta, p.tree_energy_min_frac);

%% ============ State Init ============
alive = true(p.N+1,1);
E = ones(p.N+1,1) * p.E0;

Q = cell(p.N+1,1);
for i = 1:p.N+1, Q{i} = []; end

% Metrics
gen_pkts = 0;
del_pkts = 0;
drop_pkts = 0;

sum_delay = 0;
sum_hops = 0;

energy = struct('tx',0,'rx',0,'idle',0);
diag = struct('q_overflow',0,'ttl_drop',0,'dead_drop',0,'no_route_drop',0,'prr_fail',0,'retx_exceed',0);

% Time series
time.alive_ratio = zeros(p.T,1);
time.generated = zeros(p.T,1);
time.delivered = zeros(p.T,1);

%% ============ Main Loop ============
for t = 1:p.T

    % 1) Idle energy
    E(1:p.N) = E(1:p.N) - p.Eidle;
    energy.idle = energy.idle + p.N * p.Eidle;

    alive(1:p.N) = alive(1:p.N) & (E(1:p.N) > 0);

    % 2) Random failure
    rand_fail = (rand(p.N,1) < p.p_rand);
    alive(1:p.N) = alive(1:p.N) & ~rand_fail;

    % 3) Fire
    fire_mask = is_fire_node & alive(1:p.N);
    if any(fire_mask)
        d_fire = vecnorm(pos(fire_mask,:) - p.fire_pos, 2, 2);
        p_fire = 1 - exp(-p.beta_fire .* (p.Rf - d_fire));
        p_fire = min(1, p.fire_kill_scale * p_fire);
        kill = rand(sum(fire_mask),1) < p_fire;
        idxs = find(fire_mask);
        alive(idxs(kill)) = false;
    end

    % 3.5) Periodic tree rebuild
    if startsWith(lower(p.route_mode), 'tree') && isfinite(p.tree_rebuild_period)
        if mod(t-1, max(1,round(p.tree_rebuild_period))) == 0
            Adj_tree = Adj & (alive' & alive);
            if p.USE_PRR_TH
                Adj_tree = Adj_tree & (exp_prr >= p.PRR_MIN);
            end
            parent = build_converge_tree_energy(Adj_tree, exp_prr, bs_idx, p.tree_cost, ...
                E, p.E0, p.tree_energy_beta, p.tree_energy_min_frac);
        end
    end

    % 3.6) Local parent switching (AODV simplified)
    if p.local_switch_enable
        parent = local_parent_switch(parent, alive, Adj, exp_prr, bs_idx, E, p.E0, ...
            hop_to_bs, d2bs_all, p);
    end

    % 4) Generate packets
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

    % 5) Choose senders
    nonempty = find(alive & cellfun(@(x) ~isempty(x), Q));
    nonempty(nonempty == bs_idx) = [];
    if numel(nonempty) > p.Ksend
        senders = nonempty(randperm(numel(nonempty), p.Ksend));
    else
        senders = nonempty;
    end

    % 6) Send
    for s = 1:numel(senders)
        i = senders(s);

        if isempty(Q{i}) || ~alive(i)
            continue;
        end

        pkt = Q{i}(1);
        Q{i}(1) = [];

        pkt.ttl = pkt.ttl - 1;
        if pkt.ttl <= 0
            diag.ttl_drop = diag.ttl_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        nbrs = find(Adj(i,:) & alive');
        if isempty(nbrs)
            diag.no_route_drop = diag.no_route_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        prr = exp_prr(i, nbrs) + p.noise_sigma * randn(size(nbrs));
        prr = min(1, max(p.prr_floor, prr));

        [j, prr_ij] = select_next_hop(p.route_mode, i, nbrs, prr, D(i,nbrs), ...
            bs_idx, d2bs_all, E, hop_to_bs, parent, Adj, alive);

        if j == 0
            diag.no_route_drop = diag.no_route_drop + 1;
            drop_pkts = drop_pkts + 1;
            continue;
        end

        success = false;
        for rtx = 0:p.MAX_RETX
            E(i) = E(i) - p.Etx;
            energy.tx = energy.tx + p.Etx;
            if j ~= bs_idx
                E(j) = E(j) - p.Erx;
                energy.rx = energy.rx + p.Erx;
            end

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

        pkt.hop = pkt.hop + 1;

        if j == bs_idx
            del_pkts = del_pkts + 1;
            sum_delay = sum_delay + (t - pkt.born_t);
            sum_hops = sum_hops + pkt.hop;
        else
            [Q{j}, of] = push_pkt(Q{j}, pkt, p.Qmax);
            if of
                diag.q_overflow = diag.q_overflow + 1;
                drop_pkts = drop_pkts + 1;
            end
        end
    end

    time.alive_ratio(t) = mean(alive(1:p.N));
    time.generated(t) = gen_pkts;
    time.delivered(t) = del_pkts;
end

%% ============ Summary ============
final = struct();
final.PDR = safe_div(del_pkts, gen_pkts);
final.Delay = safe_div(sum_delay, del_pkts);
final.Hops = safe_div(sum_hops, del_pkts);
final.Alive = mean(alive(1:p.N));
final.Eavg = mean(max(E(1:p.N),0));
final.Esum = sum(max(E(1:p.N),0));

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

%% ==================== Next-hop selection ====================
function [j, prr_ij] = select_next_hop(route_mode, i, nbrs, prr, d, ...
                                      bs_idx, d2bs_all, E, hop_to_bs, parent, Adj, alive)
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
        score = normalize01(E_n) + normalize01(prr);
        [~, idx] = max(score);
        j = nbrs(idx);
        prr_ij = prr(idx);

    case 'mix'
        score = 0.45*normalize01(prr) + 0.35*(1-normalize01(d2bs_n)) + 0.20*normalize01(E_n);
        [~, idx] = max(score);
        j = nbrs(idx);
        prr_ij = prr(idx);

    otherwise
        [~, idx] = min(etx);
        j = nbrs(idx);
        prr_ij = prr(idx);
end
end

%% ==================== Local parent switching ====================
function parent = local_parent_switch(parent, alive, Adj, exp_prr, bs_idx, E, E0, hop_to_bs, d2bs_all, p)
for i = 1:numel(parent)
    if i == bs_idx || ~alive(i)
        continue;
    end

    nbrs = find(Adj(i,:) & alive');
    if isempty(nbrs)
        continue;
    end

    curr = parent(i);
    curr_ok = (curr > 0) && alive(curr) && Adj(i,curr);
    need_switch = ~curr_ok;

    if curr_ok
        prr_curr = exp_prr(i, curr);
        if prr_curr < p.local_switch_prr_min
            need_switch = true;
        end
        if (E(curr) / max(E0,1e-9)) < p.local_switch_energy_min
            need_switch = true;
        end
    end

    cand = nbrs;
    if p.local_switch_hop_strict
        cand = cand(hop_to_bs(cand) < hop_to_bs(i));
    end
    if isempty(cand)
        cand = nbrs(d2bs_all(nbrs) < d2bs_all(i));
    end
    if isempty(cand)
        continue;
    end

    [best, best_cost] = select_best_parent(i, cand, exp_prr, E, E0, hop_to_bs, bs_idx, p);
    curr_cost = inf;
    if curr_ok
        curr_cost = local_parent_cost(i, curr, exp_prr, E, E0, hop_to_bs, bs_idx, p);
    end

    if need_switch || best_cost < (1 - p.local_switch_gain) * curr_cost
        parent(i) = best;
    end
end
end

function [best, best_cost] = select_best_parent(i, cand, exp_prr, E, E0, hop_to_bs, bs_idx, p)
best_cost = inf;
best = cand(1);
for k = 1:numel(cand)
    j = cand(k);
    c = local_parent_cost(i, j, exp_prr, E, E0, hop_to_bs, bs_idx, p);
    if c < best_cost
        best_cost = c;
        best = j;
    end
end
end

function cost = local_parent_cost(i, j, exp_prr, E, E0, hop_to_bs, bs_idx, p)
prr = exp_prr(i, j);
etx = 1 / max(prr, 1e-6);
if j == bs_idx
    pen = 1.0;
else
    frac = max(0, E(j)) / max(E0, 1e-9);
    frac = min(1.0, max(p.tree_energy_min_frac, frac));
    pen = 1.0 + p.tree_energy_beta * (1.0 - frac);
end
cost = etx * pen + p.local_switch_hop_weight * hop_to_bs(j);
end

%% ==================== Queue helpers ====================
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

%% ==================== Energy-aware tree ====================
function parent = build_converge_tree_energy(Adj, exp_prr, bs_idx, tree_cost, E, E0, beta, min_frac)
n = size(Adj,1);
parent = zeros(n,1);

tree_cost = lower(string(tree_cost));
use_hop = (tree_cost == "hop") || (tree_cost == "minhop");

dist = inf(n,1);
visited = false(n,1);
dist(bs_idx) = 0;

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
        base_w = 1 ./ max(p_, 1e-6);
    end

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

%% ==================== Defaults ====================
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

def.route_mode = 'tree_energy_aodv';

def.tree_cost = 'etx';
def.tree_rebuild_period = 50;
def.tree_energy_aware = true;
def.tree_energy_beta  = 2.0;
def.tree_energy_min_frac = 0.10;

def.local_switch_enable = false;
def.local_switch_prr_min = 0.20;
def.local_switch_energy_min = 0.15;
def.local_switch_gain = 0.10;
def.local_switch_hop_weight = 0.10;
def.local_switch_hop_strict = true;

def.seed = 1;
def.verbose = false;

fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(p, f) || isempty(p.(f))
        p.(f) = def.(f);
    end
end
end
