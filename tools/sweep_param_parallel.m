function res = sweep_param_parallel(route_modes, n_runs, base_seed, base, field, values, sim_func)
%SWEEP_PARAM_PARALLEL (project-aligned)
% Output matches run_topo_tree style, so plot_metric_vs_param_paper can be reused directly.
%
% res(m) contains:
%   .route_mode
%   .PDR, .PDR_std
%   .avg_delay, .avg_delay_std
%   .avg_hops,  .avg_hops_std
%   .alive_ratio, .alive_ratio_std
%   .energy_total, .energy_total_std
%   .FND, .FND_std, .HND, .HND_std, .LND, .LND_std
%
% sim_func: function handle, out = sim_func(p). If omitted, defaults to main_energy_aware_tree(p).

if nargin < 7 || isempty(sim_func)
    sim_func = @main_energy_aware_tree;
end

metrics = {'PDR','avg_delay','avg_hops','alive_ratio','energy_total','FND','HND','LND'};
res = struct();

for m = 1:numel(route_modes)
    res(m).route_mode = route_modes{m};
    for k = 1:numel(metrics)
        res(m).(metrics{k}) = zeros(numel(values),1);
        res(m).([metrics{k} '_std']) = zeros(numel(values),1);
    end
end

for m = 1:numel(route_modes)
    mode = route_modes{m};

    for vi = 1:numel(values)
        p = base;
        p.route_mode = mode;
        p.(field) = values(vi);

        % Apply energy-aware tree cfg if tool exists
        if exist('apply_energy_aware_cfg','file') == 2
            p = apply_energy_aware_cfg(p);
        end

        PDR_arr   = zeros(n_runs,1);
        Dly_arr   = zeros(n_runs,1);
        Hop_arr   = zeros(n_runs,1);
        Alive_arr = zeros(n_runs,1);
        E_arr     = zeros(n_runs,1);
        FND_arr   = zeros(n_runs,1);
        HND_arr   = zeros(n_runs,1);
        LND_arr   = zeros(n_runs,1);

        % Parallelize on runs (stable slicing)
        parfor r = 1:n_runs
            pr = p;
            pr.seed = base_seed + 100000*m + 1000*vi + r;
            rng(pr.seed, 'twister');

            out = sim_func(pr);
            f = out.final;

            PDR_arr(r) = f.PDR;

            % Map fields (support both naming styles)
            if isfield(f,'Delay')
                Dly_arr(r) = f.Delay;
            elseif isfield(f,'avg_delay')
                Dly_arr(r) = f.avg_delay;
            else
                Dly_arr(r) = NaN;
            end

            if isfield(f,'Hops')
                Hop_arr(r) = f.Hops;
            elseif isfield(f,'avg_hops')
                Hop_arr(r) = f.avg_hops;
            else
                Hop_arr(r) = NaN;
            end

            if isfield(f,'Alive')
                Alive_arr(r) = f.Alive;
            elseif isfield(f,'alive_ratio')
                Alive_arr(r) = f.alive_ratio;
            else
                Alive_arr(r) = NaN;
            end

            if isfield(f,'Esum')
                E_arr(r) = f.Esum;
            elseif isfield(f,'energy_total')
                E_arr(r) = f.energy_total;
            elseif isfield(f,'Etotal')
                E_arr(r) = f.Etotal;
            else
                E_arr(r) = NaN;
            end

            % Lifetime
            if isfield(f,'FND'), FND_arr(r) = f.FND; else, FND_arr(r) = NaN; end
            if isfield(f,'HND'), HND_arr(r) = f.HND; else, HND_arr(r) = NaN; end
            if isfield(f,'LND'), LND_arr(r) = f.LND; else, LND_arr(r) = NaN; end
        end

        % Mean + std
        res(m).PDR(vi)           = mean(PDR_arr,'omitnan');
        res(m).PDR_std(vi)       = std(PDR_arr,'omitnan');

        res(m).avg_delay(vi)     = mean(Dly_arr,'omitnan');
        res(m).avg_delay_std(vi) = std(Dly_arr,'omitnan');

        res(m).avg_hops(vi)      = mean(Hop_arr,'omitnan');
        res(m).avg_hops_std(vi)  = std(Hop_arr,'omitnan');

        res(m).alive_ratio(vi)   = mean(Alive_arr,'omitnan');
        res(m).alive_ratio_std(vi)= std(Alive_arr,'omitnan');

        res(m).energy_total(vi)  = mean(E_arr,'omitnan');
        res(m).energy_total_std(vi)= std(E_arr,'omitnan');

        res(m).FND(vi)           = mean(FND_arr,'omitnan');
        res(m).FND_std(vi)       = std(FND_arr,'omitnan');

        res(m).HND(vi)           = mean(HND_arr,'omitnan');
        res(m).HND_std(vi)       = std(HND_arr,'omitnan');

        res(m).LND(vi)           = mean(LND_arr,'omitnan');
        res(m).LND_std(vi)       = std(LND_arr,'omitnan');
    end
end
end
