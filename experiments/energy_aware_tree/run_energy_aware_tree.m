function run_energy_aware_tree()
%RUN_ENERGY_AWARE_TREE Engineering version:
% - No local functions (PARFOR-safe)
% - Traceable outputs under /results (same philosophy as run_topo_tree)
% - Adds 'tree_energy' mode while keeping the rest identical
%
% Place this file under: experiments/energy_aware_tree/
% Place helper tools under: tools/

clc; close all;

% ---------- Ensure project paths ----------
this_file = mfilename('fullpath');
this_dir  = fileparts(this_file);
proj_root = fullfile(this_dir, '..', '..');        % .../GNC
addpath(genpath(proj_root));                       % include core/tools/experiments

S = paper_plot_style();

%% ============ Base config ============
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

base.tree_cost = 'etx';
base.tree_rebuild_period = inf;

base.tree_energy_beta = 2.0;
base.tree_energy_min_frac = 0.10;

%% ============ Modes & sweeps ============
route_modes = {'tree','tree_energy','minhop','etx','greedy','energy','mindist','mix'};
topo_modes  = {'uniform','cluster','road'};
lambda_list = [0.02, 0.05, 0.08, 0.10];
n_runs = 10;

%% ============ Traceable output dir under /results ============
exp_name = 'energy_aware_tree';
if exist('make_run_dir','file') == 2
    run_dir = make_run_dir(exp_name);
else
    stamp = datestr(now,'yyyymmdd_HHMMSS');
    run_dir = fullfile(proj_root, 'results', [exp_name '_' stamp]);
    if ~exist(run_dir,'dir'), mkdir(run_dir); end
end

if exist('save_run_metadata','file') == 2
    [~, rid] = fileparts(run_dir);
    ts = datestr(now,'yyyy-mm-dd HH:MM:SS');

    meta = struct();
    meta.exp_name = exp_name;
    meta.run_id = rid;
    meta.timestamp = ts;
    meta.tag = 'energy_aware_tree';

    meta.matlab_version = version;
    meta.matlab_release = version('-release');   % <<< 关键：补齐这个字段
    meta.computer = computer;                    % <<< 关键：补齐这个字段

    % 建议补齐（可追溯更完整）
    meta.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
    meta.script = mfilename('fullpath');
    meta.note = '';

    % Git 信息（即使取不到也不影响）
    meta.git = struct('hash','','branch','','is_dirty',0);
    try
        [~, b] = system('git rev-parse --abbrev-ref HEAD'); meta.git.branch = strtrim(b);
        [~, c] = system('git rev-parse HEAD');             meta.git.hash   = strtrim(c);
        [~, d] = system('git status --porcelain');         meta.git.is_dirty = ~isempty(strtrim(d));
    catch
    end

    save_run_metadata(run_dir, meta);
end







%% ============ Parallel pool ============
if license('test','Distrib_Computing_Toolbox')
    ppool = gcp('nocreate');
    if isempty(ppool)
        parpool;
    end
end

%% ============ Run per topo_mode ============
for tm = 1:numel(topo_modes)
    topo_mode = topo_modes{tm};
    fprintf('\n================ topo_mode=%s ================\n', topo_mode);

    base_this = base;
    base_this.topo_mode = topo_mode;

    save_dir = fullfile(run_dir, topo_mode);
    if ~exist(save_dir,'dir'), mkdir(save_dir); end

    res_lambda = sweep_param_parallel(route_modes, n_runs, base_seed + 1000*tm, base_this, ...
        'lambda', lambda_list, @sim_one_run_energy_aware);

    if exist('plot_metric_vs_param_paper','file') == 2
        plot_metric_vs_param_paper(lambda_list, res_lambda.PDR_m,   route_modes, 'lambda', 'PDR',   save_dir, S);
        plot_metric_vs_param_paper(lambda_list, res_lambda.Delay_m, route_modes, 'lambda', 'Delay', save_dir, S);
        plot_metric_vs_param_paper(lambda_list, res_lambda.Hops_m,  route_modes, 'lambda', 'Hops',  save_dir, S);
        plot_metric_vs_param_paper(lambda_list, res_lambda.FND_m,   route_modes, 'lambda', 'FND',   save_dir, S);
        plot_metric_vs_param_paper(lambda_list, res_lambda.HND_m,   route_modes, 'lambda', 'HND',   save_dir, S);
        plot_metric_vs_param_paper(lambda_list, res_lambda.LND_m,   route_modes, 'lambda', 'LND',   save_dir, S);
    end

    save(fullfile(save_dir, 'results_energy_aware_tree.mat'), ...
        'route_modes','lambda_list','topo_mode','n_runs','base_this','res_lambda');
end

fprintf('\nAll done. Results saved to: %s\n', run_dir);
end
