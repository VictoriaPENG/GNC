function run_baseline_experiments_parfor()
% ============================================================
% run_baseline_experiments.m  (Parallel version using parfor)
% ------------------------------------------------------------
% 一键生成"基线对比曲线"（论文规范版 + 并行加速）
%
% 依赖文件（同目录）：
%   - main_sim.m
%   - paper_plot_style.m
%   - plot_metric_vs_param_paper.m
%
% 输出：
%   ./figs/ 目录下生成多张图（pdf/eps/png）
%   ./figs/results_baseline.mat 保存所有结果结构体
% ============================================================

close all; clc;

%% ============ 多次重复（用于均值+95%CI） ============
n_runs = 20;          % 有并行后建议 8~20
base_seed = 1;

%% ============ 输出目录 ============
outdir = sprintf('./figs_n_run_%d', n_runs);
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ============ 路由策略集合（baseline） ============
route_modes = {'minhop','etx','greedy','energy','mindist','mix'};



%% ============ 并行池（自动开启） ============
pool = gcp('nocreate');
if isempty(pool)
    parpool;  %#ok<*NOPRT>
end

%% ============ 基础参数（固定项） ============
base = struct();
base.verbose = true; %仿真日志输出控制
base.seed = base_seed;

% 网络与仿真
base.Lx = 1000; base.Ly = 1000;
base.N  = 700;
base.Rc = 150;
base.T  = 2000;

% 基站/火区
base.BS_pos   = [500, 500];
base.fire_pos = [250, 250];
base.Rf       = 220;

% 业务负载（会做 sweep）
base.lambda = 0.05;

% 信道（会做 sweep）
base.alpha = 0.03;
base.prr_floor = 0.05;
base.noise_sigma = 0.02;

% 链路阈值（可选）
base.USE_PRR_TH = true;
base.PRR_MIN = 0.15;

% 能耗模型
base.E0    = 100;
base.Etx   = 0.9;
base.Erx   = 0.4;
base.Eidle = 0.001;

% 失效模型
base.p_rand = 1e-5;
base.beta_fire = 0.05;
base.fire_kill_scale = 1e-4;

% 队列/转发/MAC
base.Qmax = 50;
base.TTLmax = 25;
base.Ksend = 20;
base.MAX_RETX = 5;

% 重传放队尾（推荐）
base.RETX_TO_TAIL = true;

%% ============ 论文级绘图风格 ============
S = paper_plot_style();

%% ============ sweep 配置 ============
lambda_list = [0.01 0.02 0.05 0.08 0.10];
alpha_list  = [0.005 0.01 0.02 0.03 0.04];
N_list      = [200 400 700 1000];
Rc_list     = [100 120 150 180 220];

%% ============ 开始跑实验（并行） ============
fprintf("Running baselines (parallel) with n_runs=%d ...\n", n_runs);

res_lambda = sweep_param_parallel(route_modes, n_runs, base_seed, base, 'lambda', lambda_list);
res_alpha  = sweep_param_parallel(route_modes, n_runs, base_seed, base, 'alpha',  alpha_list);
res_N      = sweep_param_parallel(route_modes, n_runs, base_seed, base, 'N',      N_list);
res_Rc     = sweep_param_parallel(route_modes, n_runs, base_seed, base, 'Rc',     Rc_list);

%% ============ 画论文级曲线（均值+95%CI） ============
plot_metric_vs_param_paper(res_lambda, lambda_list, ...
    'PDR', 'PDR', '\lambda (packet generation probability)', ...
    fullfile(outdir,'PDR_vs_lambda'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_lambda, lambda_list, ...
    'avg_delay', 'Average Delay (steps)', '\lambda (packet generation probability)', ...
    fullfile(outdir,'Delay_vs_lambda'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_alpha, alpha_list, ...
    'PDR', 'PDR', '\alpha (path-loss factor)', ...
    fullfile(outdir,'PDR_vs_alpha'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_N, N_list, ...
    'avg_delay', 'Average Delay (steps)', 'N (number of sensors)', ...
    fullfile(outdir,'Delay_vs_N'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_N, N_list, ...
    'energy_total', 'Total Energy (idle+tx+rx)', 'N (number of sensors)', ...
    fullfile(outdir,'Energy_vs_N'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_Rc, Rc_list, ...
    'PDR', 'PDR', 'R_c (communication radius, m)', ...
    fullfile(outdir,'PDR_vs_Rc'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

%% ============ AliveRatio vs Time（代表性单次） ============
figure('Color','w','Units','inches', 'Position',[1 1 S.figureW_in S.figureH_in]);
ax = axes(); hold(ax,'on'); grid(ax,'on');
set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, 'LineWidth', 0.8, 'Box','on');
ax.GridAlpha = S.gridAlpha;

t = (1:base.T)';
markers = {'o','s','^','d','v','>','<','p','h','x','+'};

for m = 1:numel(route_modes)
    p = base;
    p.route_mode = route_modes{m};
    p.seed = base_seed + 1000*m;
    rng(p.seed, 'twister'); % 保证这次代表性运行可复现
    out = main_sim_test(p);

    mk = markers{1 + mod(m-1, numel(markers))};
    plot(ax, t, out.time.alive_ratio, 'LineWidth', S.lineWidth, ...
        'Marker', mk, 'MarkerSize', S.markerSize, ...
        'DisplayName', route_modes{m});
end

xlabel(ax, 'Time step');
ylabel(ax, 'Alive Ratio');
legend(ax, route_modes, 'Location', S.legendLocation, 'Box', S.legendBox, ...
    'FontName', S.fontName, 'FontSize', S.fontSize-1);

tight_inset(ax);
export_figure(gcf, fullfile(outdir,'AliveRatio_vs_Time'), S);

%% ============ 寿命指标（FND/HND/LND）并行统计 + 条形图 ============
life = sweep_single_point_parallel(route_modes, n_runs, base_seed, base);
plot_lifetime_bar(life, n_runs, S, fullfile(outdir,'Lifetime_FND_HND_LND'));

%% ============ 保存所有结果 ============
results = struct();
results.base = base;
results.route_modes = route_modes;
results.n_runs = n_runs;
results.res_lambda = res_lambda;
results.res_alpha  = res_alpha;
results.res_N      = res_N;
results.res_Rc     = res_Rc;
results.life       = life;

save(fullfile(outdir,'results_baseline.mat'), '-struct', 'results');
fprintf("\nDone. Figures and results saved to: %s\n", outdir);

end

%% =======================================================================
% Parallel sweep: 对单个参数做 sweep，并输出均值/标准差（用于CI绘图）
% 并行粒度：对 n_runs 做 parfor（最稳定，避免嵌套 parfor）
% =======================================================================
function res = sweep_param_parallel(route_modes, n_runs, base_seed, base, field, values)

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
        v = values(vi);

        % --- 用普通数组承接 parfor 输出（parfor 可切片） ---
        PDR_arr         = zeros(n_runs,1);
        delay_arr       = zeros(n_runs,1);
        hops_arr        = zeros(n_runs,1);
        alive_arr       = zeros(n_runs,1);
        energy_arr      = zeros(n_runs,1);
        FND_arr         = zeros(n_runs,1);
        HND_arr         = zeros(n_runs,1);
        LND_arr         = zeros(n_runs,1);

        parfor r = 1:n_runs
            p = base;
            p.route_mode = mode;
            p.(field) = v;

            % 每次 run 的确定性 seed（可复现）
            p.seed = base_seed + 10000*m + 100*vi + r;
            rng(p.seed, 'twister');

            out = main_sim_test(p);
            f = out.final;

            PDR_arr(r)    = f.PDR;
            delay_arr(r)  = f.avg_delay;
            hops_arr(r)   = f.avg_hops;
            alive_arr(r)  = f.alive_ratio;
            energy_arr(r) = f.energy_total;

            FND_arr(r) = nan_to_T(f.FND, p.T);
            HND_arr(r) = nan_to_T(f.HND, p.T);
            LND_arr(r) = nan_to_T(f.LND, p.T);
        end

        % --- 汇总到 res（均值/标准差） ---
        res(m).PDR(vi)           = mean(PDR_arr);
        res(m).PDR_std(vi)       = std(PDR_arr);

        res(m).avg_delay(vi)     = mean(delay_arr);
        res(m).avg_delay_std(vi) = std(delay_arr);

        res(m).avg_hops(vi)      = mean(hops_arr);
        res(m).avg_hops_std(vi)  = std(hops_arr);

        res(m).alive_ratio(vi)       = mean(alive_arr);
        res(m).alive_ratio_std(vi)   = std(alive_arr);

        res(m).energy_total(vi)      = mean(energy_arr);
        res(m).energy_total_std(vi)  = std(energy_arr);

        res(m).FND(vi)           = mean(FND_arr);
        res(m).FND_std(vi)       = std(FND_arr);

        res(m).HND(vi)           = mean(HND_arr);
        res(m).HND_std(vi)       = std(HND_arr);

        res(m).LND(vi)           = mean(LND_arr);
        res(m).LND_std(vi)       = std(LND_arr);
    end
end

end


%% =======================================================================
% 单点重复：寿命统计（并行）
% =======================================================================
function life = sweep_single_point_parallel(route_modes, n_runs, base_seed, base)

life = struct();
for m = 1:numel(route_modes)
    mode = route_modes{m};

    FND = zeros(n_runs,1);
    HND = zeros(n_runs,1);
    LND = zeros(n_runs,1);

    parfor r = 1:n_runs
        p = base;
        p.route_mode = mode;
        p.seed = base_seed + 55500*m + r;
        rng(p.seed, 'twister');

        out = main_sim_test(p);
        f = out.final;

        FND(r) = nan_to_T(f.FND, p.T);
        HND(r) = nan_to_T(f.HND, p.T);
        LND(r) = nan_to_T(f.LND, p.T);
    end

    life(m).route_mode = mode;
    life(m).FND_mean = mean(FND); life(m).FND_std = std(FND);
    life(m).HND_mean = mean(HND); life(m).HND_std = std(HND);
    life(m).LND_mean = mean(LND); life(m).LND_std = std(LND);
end

end

%% =======================================================================
% 寿命条形图（FND/HND/LND）
% =======================================================================
function plot_lifetime_bar(life, n_runs, S, savebase)

modes = {life.route_mode};

FND = [life.FND_mean]'; FNDs = [life.FND_std]';
HND = [life.HND_mean]'; HNDs = [life.HND_std]';
LND = [life.LND_mean]'; LNDs = [life.LND_std]';

z = 1.96;
FNDci = z * (FNDs ./ sqrt(n_runs));
HNDci = z * (HNDs ./ sqrt(n_runs));
LNDci = z * (LNDs ./ sqrt(n_runs));

Y = [FND HND LND];
CI = [FNDci HNDci LNDci];

figure('Color','w','Units','inches', 'Position',[1 1 max(5.5,S.figureW_in*1.6) S.figureH_in]);
ax = axes(); hold(ax,'on'); grid(ax,'on');
set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, 'LineWidth', 0.8, 'Box','on');
ax.GridAlpha = S.gridAlpha;

hb = bar(ax, Y); 
ng = size(Y,1);
nb = size(Y,2);

xpos = nan(ng, nb);
for b = 1:nb
    xpos(:,b) = hb(b).XEndPoints;
end
for b = 1:nb
    errorbar(ax, xpos(:,b), Y(:,b), CI(:,b), 'k', 'LineStyle','none', 'CapSize',6);
end

set(ax, 'XTick', 1:ng, 'XTickLabel', modes);
ylabel(ax, 'Time step');
legend(ax, {'FND','HND','LND'}, 'Location', S.legendLocation, 'Box', S.legendBox, ...
    'FontName', S.fontName, 'FontSize', S.fontSize-1);

tight_inset(ax);
export_figure(gcf, savebase, S);

end

%% =======================================================================
% 小工具
% =======================================================================
function y = nan_to_T(x, T)
if isnan(x), y = T; else, y = x; end
end

function tight_inset(ax)
outer = ax.OuterPosition;
ti = ax.TightInset;
left = outer(1) + ti(1);
bottom = outer(2) + ti(2);
ax_width = outer(3) - ti(1) - ti(3);
ax_height = outer(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end

function export_figure(fig, savebase, S)
set(fig, 'PaperPositionMode','auto');
print(fig, [savebase '.pdf'], '-dpdf', '-vector');
print(fig, [savebase '.eps'], '-depsc', '-vector');
print(fig, [savebase '.png'], ['-r' num2str(S.exportDPI)], '-dpng');
end
