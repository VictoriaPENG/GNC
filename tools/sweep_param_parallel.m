function res = sweep_param_parallel(route_modes, n_runs, base_seed, base_cfg, param_name, param_list, sim_func)
%SWEEP_PARAM_PARALLEL Parallel sweep utility (PARFOR-safe, engineered).
%
% res = sweep_param_parallel(route_modes, n_runs, base_seed, base_cfg, ...
%                            param_name, param_list, sim_func)
%
% Output arrays are [n_mode x n_param x n_runs]

if nargin < 7 || isempty(sim_func)
    error('sweep_param_parallel:MissingSimFunc', 'sim_func must be provided.');
end

n_mode  = numel(route_modes);
n_param = numel(param_list);
total   = n_mode * n_param * n_runs;

% ---- Use flat arrays for PARFOR slicing safety ----
PDR_flat   = nan(total,1);
Delay_flat = nan(total,1);
Hops_flat  = nan(total,1);
FND_flat   = nan(total,1);
HND_flat   = nan(total,1);
LND_flat   = nan(total,1);

parfor idx = 1:total
    [m, pi, r] = ind2sub([n_mode, n_param, n_runs], idx);

    p = base_cfg;
    p.route_mode = route_modes{m};
    p.(param_name) = param_list(pi);

    % Deterministic seed per (mode, param, run)
    p.seed = base_seed + 10000*m + 100*pi + r;
    rng(p.seed, 'twister');

    % Apply energy-aware config if tool exists (no hard dependency)
    if exist('apply_energy_aware_cfg','file') == 2
        p = apply_energy_aware_cfg(p);
    end

    out = sim_func(p);
    f = out.final;

    PDR_flat(idx)   = f.PDR;
    Delay_flat(idx) = f.Delay;
    Hops_flat(idx)  = f.Hops;
    FND_flat(idx)   = f.FND;
    HND_flat(idx)   = f.HND;
    LND_flat(idx)   = f.LND;
end

% ---- Reshape back to 3D ----
PDR   = reshape(PDR_flat,   [n_mode, n_param, n_runs]);
Delay = reshape(Delay_flat, [n_mode, n_param, n_runs]);
Hops  = reshape(Hops_flat,  [n_mode, n_param, n_runs]);
FND   = reshape(FND_flat,   [n_mode, n_param, n_runs]);
HND   = reshape(HND_flat,   [n_mode, n_param, n_runs]);
LND   = reshape(LND_flat,   [n_mode, n_param, n_runs]);

res = struct();
res.PDR   = PDR;
res.Delay = Delay;
res.Hops  = Hops;
res.FND   = FND;
res.HND   = HND;
res.LND   = LND;

% Means
res.PDR_m   = mean(PDR,   3, 'omitnan');
res.Delay_m = mean(Delay, 3, 'omitnan');
res.Hops_m  = mean(Hops,  3, 'omitnan');
res.FND_m   = mean(FND,   3, 'omitnan');
res.HND_m   = mean(HND,   3, 'omitnan');
res.LND_m   = mean(LND,   3, 'omitnan');
end
