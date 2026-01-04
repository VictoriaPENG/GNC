function run_energy_aware_tree()
% ============================================================
% run_topo_tree.m
% ------------------------------------------------------------
% 基线对比实验（论文版）+ 支持拓扑模式（uniform/cluster/road）+ 汇聚树
% 调用 main_energy_aware_tree.m 进行仿真
%
% [新增] tree_energy：
%   - 能量感知建树（Energy-aware Tree）
%   - 在 run 中自动开启 p.tree_energy_aware=true，并设置重建周期
% ============================================================

clc; close all;

%% ============ 论文级绘图风格 ============
S = paper_plot_style();

%% ============ 基础参数（固定项） ============
base_seed = 1;
base = struct();
base.verbose = true;
base.seed = base_seed;

% 网络与仿真
base.Lx = 1000; base.Ly = 1000;
base.N  = 700;
base.Rc = 150;
base.T  = 2000;

% 基站/火区
base.BS_pos   = [500, 500];
base.fire_pos = [250, 250];
base.Rf = 220;

% 信道
base.alpha = 0.01;
base.prr_floor = 0.08;
base.noise_sigma = 0.02;

% PRR 阈值
base.USE_PRR_TH = true;
base.PRR_MIN = 0.15;

% 能耗
base.E0 = 100;
base.Etx = 0.9;
base.Erx = 0.4;
base.Eidle = 0.001;

% 失效/火灾
base.p_rand = 1e-5;
base.beta_fire = 0.05;
base.fire_kill_scale = 1e-4;

% MAC/队列
base.Qmax = 50;
base.TTLmax = 25;
base.Ksend = 20;
base.MAX_RETX = 5;

% 汇聚树参数（默认）
base.tree_cost = 'etx';
base.tree_rebuild_period = inf;

% 能量感知建树参数（可调）
base.tree_energy_beta = 2.0;
base.tree_energy_min_frac = 0.10;

%% ============ 路由模式（baseline） ============
route_modes = {'tree','tree_energy','minhop','etx','greedy','energy','mindist','mix'};

%% ============ sweep：topo_mode ============
topo_modes = {'uniform','cluster','road'};

% sweep：lambda（示例）
lambda_list = [0.02, 0.05, 0.08, 0.10];

n_run = 20;   % 你需要更稳就调大
out_dir = fullfile(pwd, 'figs_topo_article');

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

for tm = 1:numel(topo_modes)
    topo_mode = topo_modes{tm};

    fprintf('\n================ topo_mode=%s ================\n', topo_mode);

    save_dir = fullfile(out_dir, topo_mode);
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    % 结果数组：mode x lambda x run
    PDR_arr = nan(numel(route_modes), numel(lambda_list), n_run);
    Delay_arr = nan(size(PDR_arr));
    Hops_arr = nan(size(PDR_arr));
    FND_arr = nan(size(PDR_arr));
    HND_arr = nan(size(PDR_arr));
    LND_arr = nan(size(PDR_arr));

    for m = 1:numel(route_modes)
        mode = route_modes{m};

        for li = 1:numel(lambda_list)
            lam = lambda_list(li);

            for r = 1:n_run
                p = base;
                p.topo_mode = topo_mode;
                p.lambda = lam;

                p.route_mode = mode;

                % 能量感知建树：仅对 tree_energy 模式启用
                if strcmpi(p.route_mode, 'tree_energy')
                    p.tree_energy_aware = true;
                    p.tree_rebuild_period = 50; % 建议周期性重建以跟踪能量变化
                else
                    p.tree_energy_aware = false;
                end

                % 固定随机性（确保可复现：不同 run 改 seed）
                p.seed = base_seed + 1000*tm + 100*m + r;

                rng(p.seed, 'twister');

                out = main_energy_aware_tree(p);
                f = out.final;

                PDR_arr(m, li, r)   = f.PDR;
                Delay_arr(m, li, r) = f.Delay;
                Hops_arr(m, li, r)  = f.Hops;
                FND_arr(m, li, r)   = f.FND;
                HND_arr(m, li, r)   = f.HND;
                LND_arr(m, li, r)   = f.LND;
            end
        end
    end

    % 取均值
    PDR_m = mean(PDR_arr, 3, 'omitnan');
    Delay_m = mean(Delay_arr, 3, 'omitnan');
    Hops_m = mean(Hops_arr, 3, 'omitnan');
    FND_m = mean(FND_arr, 3, 'omitnan');
    HND_m = mean(HND_arr, 3, 'omitnan');
    LND_m = mean(LND_arr, 3, 'omitnan');

    % 画图
    plot_metric_vs_param_paper(lambda_list, PDR_m, route_modes, 'lambda', 'PDR', save_dir, S);
    plot_metric_vs_param_paper(lambda_list, Delay_m, route_modes, 'lambda', 'Delay', save_dir, S);
    plot_metric_vs_param_paper(lambda_list, Hops_m, route_modes, 'lambda', 'Hops', save_dir, S);
    plot_metric_vs_param_paper(lambda_list, FND_m, route_modes, 'lambda', 'FND', save_dir, S);
    plot_metric_vs_param_paper(lambda_list, HND_m, route_modes, 'lambda', 'HND', save_dir, S);
    plot_metric_vs_param_paper(lambda_list, LND_m, route_modes, 'lambda', 'LND', save_dir, S);

    % 保存 mat
    save(fullfile(save_dir, 'results_tree_energy.mat'), ...
        'route_modes','lambda_list','topo_mode', ...
        'PDR_arr','Delay_arr','Hops_arr','FND_arr','HND_arr','LND_arr', ...
        'PDR_m','Delay_m','Hops_m','FND_m','HND_m','LND_m');
end

fprintf('\nAll done. Results saved to: %s\n', out_dir);
end
