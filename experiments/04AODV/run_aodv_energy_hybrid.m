function run_aodv_energy_hybrid()
% ============================================================
% Hybrid routing experiment:
%   - Dynamic energy-aware tree (periodic rebuild)
%   - AODV simplified local parent switching embedded
% ------------------------------------------------------------
% Output: results/<run_id>/aodv_energy_hybrid/<topo_mode>/...
% ============================================================

clc; close all;

% Tag from env
 tag = getenv('GNC_TAG');
if isempty(tag), tag = 'aodv_energy_hybrid'; end

run_aodv_energy_hybrid_impl(tag);

end

function run_aodv_energy_hybrid_impl(tag)
S = paper_plot_style();

route_modes = {'tree_energy_dyn50', 'tree_energy_aodv_dyn50'};
topo_modes = {'uniform','cluster','road'};

n_runs = 20;
base_seed = 1;

base = struct();
base.verbose = true;
base.seed = base_seed;

base.Lx = 1000; base.Ly = 1000;
base.N  = 700;
base.Rc = 150;
base.T  = 2000;

base.BS_pos   = [500, 500];
base.fire_pos = [250, 250];
base.Rf = 220;

base.lambda = 0.02;
base.alpha = 0.01;
base.prr_floor = 0.08;
base.noise_sigma = 0.02;

base.USE_PRR_TH = true;
base.PRR_MIN = 0.15;

base.E0 = 100;
base.Etx = 0.9;
base.Erx = 0.4;
base.Eidle = 0.001;

base.p_rand = 1e-5;
base.beta_fire = 0.05;
base.fire_kill_scale = 1e-4;

base.Qmax = 50;
base.TTLmax = 25;
base.Ksend = 20;
base.MAX_RETX = 5;
base.RETX_TO_TAIL = true;

base.hotspot_enable = false;
base.hotspot_ratio  = 0.30;
base.hotspot_sigma  = 40;

base.Kc = 8;
base.cluster_sigma = 60;
base.cluster_center_mode = 'random';
base.cluster_mix_uniform_ratio = 0.15;

base.road_count = 3;
base.road_sigma = 18;
base.road_lane_width = 18;
base.road_uniform_mix_ratio = 0.20;

base.tree_cost = 'etx';
base.tree_rebuild_period = 50;
base.tree_energy_beta = 2.0;
base.tree_energy_min_frac = 0.10;

base.local_switch_prr_min = 0.20;
base.local_switch_energy_min = 0.15;
base.local_switch_gain = 0.10;
base.local_switch_hop_weight = 0.10;
base.local_switch_hop_strict = true;

lambda_list = [0.01 0.02 0.05 0.08 0.10];

[run_dir, meta] = make_run_dir('aodv_energy_hybrid', 'Tag', tag, 'Script', mfilename('fullpath'));
root_outdir = fullfile(run_dir, 'aodv_energy_hybrid');
if ~exist(root_outdir,'dir'), mkdir(root_outdir); end

cfg = struct();
cfg.route_modes = route_modes;
cfg.n_runs = n_runs;
cfg.base_seed = base_seed;
cfg.base = base;
cfg.topo_modes = topo_modes;
cfg.lambda_list = lambda_list;

for tm = 1:numel(topo_modes)
    topo_mode = topo_modes{tm};

    outdir = fullfile(root_outdir, topo_mode);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    base_this = base;
    base_this.topo_mode = topo_mode;

    res_lambda = sweep_param_parallel_hybrid(route_modes, n_runs, base_seed, base_this, 'lambda', lambda_list);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'PDR', 'PDR', '\lambda (packet generation probability)', ...
        fullfile(outdir,'PDR_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'avg_delay', 'Average delay', '\lambda (packet generation probability)', ...
        fullfile(outdir,'Delay_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'avg_hops', 'Average hops', '\lambda (packet generation probability)', ...
        fullfile(outdir,'Hops_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'FND', 'FND', '\lambda (packet generation probability)', ...
        fullfile(outdir,'FND_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'HND', 'HND', '\lambda (packet generation probability)', ...
        fullfile(outdir,'HND_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    plot_metric_vs_param_paper(res_lambda, lambda_list, ...
        'LND', 'LND', '\lambda (packet generation probability)', ...
        fullfile(outdir,'LND_vs_lambda'), n_runs, 'Style', S, 'CI', true, 'Legend', true);

    save(fullfile(outdir,'results_aodv_energy_hybrid.mat'), ...
        'route_modes','topo_mode','n_runs','base_this','lambda_list','res_lambda');
end

summary = struct();
summary.topo_modes = topo_modes;
summary.route_modes = route_modes;
summary.n_runs = n_runs;

save_run_metadata(root_outdir, meta, 'Config', cfg, 'ResultsSummary', summary);
fprintf("ALL DONE. Hybrid outputs saved under: %s\n", root_outdir);

end

function res = sweep_param_parallel_hybrid(route_modes, n_runs, base_seed, base, field, values)
res = struct();
for m = 1:numel(route_modes)
    res(m).route_mode = route_modes{m};
    res(m).PDR = zeros(numel(values),1);
    res(m).PDR_std = zeros(numel(values),1);
    res(m).avg_delay = zeros(numel(values),1);
    res(m).avg_delay_std = zeros(numel(values),1);
    res(m).avg_hops = zeros(numel(values),1);
    res(m).avg_hops_std = zeros(numel(values),1);
    res(m).alive_ratio = zeros(numel(values),1);
    res(m).alive_ratio_std = zeros(numel(values),1);
    res(m).energy_total = zeros(numel(values),1);
    res(m).energy_total_std = zeros(numel(values),1);
    res(m).FND = zeros(numel(values),1); res(m).FND_std = zeros(numel(values),1);
    res(m).HND = zeros(numel(values),1); res(m).HND_std = zeros(numel(values),1);
    res(m).LND = zeros(numel(values),1); res(m).LND_std = zeros(numel(values),1);
end

for m = 1:numel(route_modes)
    mode = route_modes{m};

    for vi = 1:numel(values)
        p = base;
        p.(field) = values(vi);

        [p, route_mode_name] = map_mode_to_config(p, mode);
        p.route_mode = route_mode_name;

        PDR_arr = zeros(n_runs,1);
        delay_arr = zeros(n_runs,1);
        hops_arr = zeros(n_runs,1);
        alive_arr = zeros(n_runs,1);
        energy_arr = zeros(n_runs,1);
        FND_arr = zeros(n_runs,1);
        HND_arr = zeros(n_runs,1);
        LND_arr = zeros(n_runs,1);

        parfor r = 1:n_runs
            pr = p;
            pr.seed = base_seed + 10000*m + 100*vi + r;
            rng(pr.seed,'twister');

            out = main_aodv_energy_hybrid(pr);
            f = out.final;

            PDR_arr(r) = f.PDR;

            if isfield(f,'avg_delay'), delay_arr(r)=f.avg_delay;
            elseif isfield(f,'Delay'), delay_arr(r)=f.Delay;
            else, delay_arr(r)=NaN; end

            if isfield(f,'avg_hops'), hops_arr(r)=f.avg_hops;
            elseif isfield(f,'Hops'), hops_arr(r)=f.Hops;
            else, hops_arr(r)=NaN; end

            if isfield(f,'alive_ratio'), alive_arr(r)=f.alive_ratio;
            elseif isfield(f,'Alive'), alive_arr(r)=f.Alive;
            else, alive_arr(r)=NaN; end

            if isfield(f,'energy_total'), energy_arr(r)=f.energy_total;
            elseif isfield(f,'Etotal'), energy_arr(r)=f.Etotal;
            elseif isfield(f,'Esum'), energy_arr(r)=f.Esum;
            else, energy_arr(r)=NaN; end

            FND_arr(r) = nan_to_T(f.FND, pr.T);
            HND_arr(r) = nan_to_T(f.HND, pr.T);
            LND_arr(r) = nan_to_T(f.LND, pr.T);
        end

        res(m).PDR(vi) = mean(PDR_arr); res(m).PDR_std(vi) = std(PDR_arr);
        res(m).avg_delay(vi) = mean(delay_arr); res(m).avg_delay_std(vi) = std(delay_arr);
        res(m).avg_hops(vi) = mean(hops_arr); res(m).avg_hops_std(vi) = std(hops_arr);
        res(m).alive_ratio(vi) = mean(alive_arr); res(m).alive_ratio_std(vi) = std(alive_arr);
        res(m).energy_total(vi) = mean(energy_arr); res(m).energy_total_std(vi) = std(energy_arr);

        res(m).FND(vi)=mean(FND_arr); res(m).FND_std(vi)=std(FND_arr);
        res(m).HND(vi)=mean(HND_arr); res(m).HND_std(vi)=std(HND_arr);
        res(m).LND(vi)=mean(LND_arr); res(m).LND_std(vi)=std(LND_arr);
    end
end
end

function [p, route_mode_name] = map_mode_to_config(p, mode)
mode = lower(string(mode));
route_mode_name = 'tree_energy';

p.tree_energy_aware = true;

period = parse_rebuild_period(mode);
if isfinite(period)
    p.tree_rebuild_period = period;
end

if contains(mode, "aodv")
    route_mode_name = 'tree_energy_aodv';
    p.local_switch_enable = true;
else
    route_mode_name = 'tree_energy';
    p.local_switch_enable = false;
end
end

function period = parse_rebuild_period(mode)
period = NaN;
expr = regexp(mode, "dyn(\\d+)", "tokens", "once");
if ~isempty(expr)
    period = str2double(expr{1});
end
if ~isfinite(period)
    period = 50;
end
end

function y = nan_to_T(x, T)
if isnan(x)
    y = T;
else
    y = x;
end
end
