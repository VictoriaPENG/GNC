function run_baseline_for_article()
% ============================================================
% run_baseline_experiments.m
% ------------------------------------------------------------
% 一键生成"基线对比曲线"（论文规范版）
%
% 依赖文件（同目录）：
%   - main_sim.m                    % 你的仿真主函数（返回 out.final/out.time）
%   - paper_plot_style.m            % 统一论文绘图风格
%   - plot_metric_vs_param_paper.m  % 论文级绘图（均值+95%CI阴影带，导出pdf/eps/png）
%
% 输出：
%   ./figs/ 目录下生成多张图（pdf/eps/png）
%   ./figs/results_baseline.mat 保存所有结果结构体
% ============================================================

close all; clc;

%% ============ 输出目录 ============
outdir = 'figs';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ============ 路由策略集合（baseline） ============
route_modes = {'minhop','etx','greedy','energy','mindist','mix'};

%% ============ 多次重复（用于均值+95%CI） ============
n_runs = 5;          % 建议 5~10
base_seed = 1;       % 改这个可整体换随机序列

%% ============ 基础参数（固定项） ============
% 这些参数会被 sweep 覆盖：lambda/alpha/N/Rc 等
base = struct();
base.verbose = false;       % 批量实验不要刷屏
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

% 建议：先关阈值做纯对比；如果要提升可用性可打开
base.USE_PRR_TH = false;
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

% 推荐：重传放队尾，避免差链路包"卡死"
base.RETX_TO_TAIL = true;

%% ============ 论文级绘图风格 ============
S = paper_plot_style();

% 可追溯输出目录（Run ID + Tag）
tag = getenv('GNC_TAG');
if isempty(tag), tag = 'baseline'; end

exp_name = 'baseline';
[run_dir, meta] = make_run_dir(exp_name, 'Tag', tag, 'Script', mfilename('fullpath'));
outdir = run_dir;

cfg = struct();
cfg.route_modes = route_modes;
cfg.n_runs = n_runs;
cfg.base_seed = base_seed;
cfg.base = base;
cfg.sweeps = struct('lambda', lambda_list, 'alpha', alpha_list, 'N', N_list, 'Rc', Rc_list);
save_run_metadata(outdir, meta, 'Config', cfg);

%% ============ sweep 配置 ============
% 1) PDR/Delay vs lambda（负载鲁棒性）
lambda_list = [0.01 0.02 0.05 0.08 0.10];

% 2) PDR vs alpha（信道恶化鲁棒性）
alpha_list = [0.005 0.01 0.02 0.03 0.04];

% 3) Delay/Energy vs N（规模性）
N_list = [200 400 700 1000];

% 4) 可选：PDR vs Rc（网络稠密度）
Rc_list = [100 120 150 180 220];

%% ============ 开始跑实验 ============
fprintf("Running baselines with n_runs=%d ...\n", n_runs);

res_lambda = sweep_param(route_modes, n_runs, base_seed, base, 'lambda', lambda_list);
res_alpha  = sweep_param(route_modes, n_runs, base_seed, base, 'alpha',  alpha_list);
res_N      = sweep_param(route_modes, n_runs, base_seed, base, 'N',      N_list);
res_Rc     = sweep_param(route_modes, n_runs, base_seed, base, 'Rc',     Rc_list);

%% ============ 画论文级曲线（均值+95%CI） ============
% --- PDR/Delay vs lambda ---
plot_metric_vs_param_paper(res_lambda, lambda_list, ...
    'PDR', 'PDR', '\lambda (packet generation probability)', ...
    fullfile(outdir,'PDR_vs_lambda'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_lambda, lambda_list, ...
    'avg_delay', 'Average Delay (steps)', '\lambda (packet generation probability)', ...
    fullfile(outdir,'Delay_vs_lambda'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

% --- PDR vs alpha ---
plot_metric_vs_param_paper(res_alpha, alpha_list, ...
    'PDR', 'PDR', '\alpha (path-loss factor)', ...
    fullfile(outdir,'PDR_vs_alpha'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

% --- Delay/Energy vs N ---
plot_metric_vs_param_paper(res_N, N_list, ...
    'avg_delay', 'Average Delay (steps)', 'N (number of sensors)', ...
    fullfile(outdir,'Delay_vs_N'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

plot_metric_vs_param_paper(res_N, N_list, ...
    'energy_total', 'Total Energy (idle+tx+rx)', 'N (number of sensors)', ...
    fullfile(outdir,'Energy_vs_N'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

% --- PDR vs Rc ---
plot_metric_vs_param_paper(res_Rc, Rc_list, ...
    'PDR', 'PDR', 'R_c (communication radius, m)', ...
    fullfile(outdir,'PDR_vs_Rc'), n_runs, ...
    'Style', S, 'CI', true, 'Legend', true);

%% ============ AliveRatio vs Time（寿命曲线：代表性单次） ============
% 注：寿命曲线如果也要均值阴影带，需要保存每次 run 的时间序列（成本更高）
% 这里先给"代表性 run"用于论文 Fig.（足够当 baseline 结果展示）
figure('Color','w','Units','inches', 'Position',[1 1 S.figureW_in S.figureH_in]);
ax = axes(); hold(ax,'on'); grid(ax,'on');
set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, 'LineWidth', 0.8, 'Box','on');
ax.GridAlpha = S.gridAlpha;

t = (1:base.T)';
markers = {'o','s','^','d','v','>','<','p','h','x','+'};

for m = 1:numel(route_modes)
    p = base;
    p.route_mode = route_modes{m};
    p.seed = base_seed + 1000*m;  % 每条曲线不同 seed
    out = main_baseline(p);

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

%% ============ 导出寿命指标对比（FND/HND/LND vs route） ============
% 用 lambda=base.lambda、alpha=base.alpha、N=base.N 的默认工况：
life = sweep_single_point(route_modes, n_runs, base_seed, base);

% 画条形图（均值±95%CI）
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
results.meta       = meta;
results.cfg        = cfg;

save(fullfile(outdir,'results_baseline.mat'), '-struct', 'results');
fprintf("\nDone. Figures and results saved to: %s\n", outdir);

end

%% =======================================================================
% Sweep: 对单个参数做 sweep，并输出均值/标准差（用于CI绘图）
% =======================================================================
function res = sweep_param(route_modes, n_runs, base_seed, base, field, values)
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

        tmp = struct();
        for k = 1:numel(metrics)
            tmp.(metrics{k}) = zeros(n_runs,1);
        end

        for r = 1:n_runs
            p = base;
            p.route_mode = mode;
            p.(field) = v;
            p.seed = base_seed + 10000*m + 100*vi + r;

            out = main_baseline(p);
            f = out.final;

            tmp.PDR(r)        = f.PDR;
            tmp.avg_delay(r)  = f.avg_delay;
            tmp.avg_hops(r)   = f.avg_hops;
            tmp.alive_ratio(r)= f.alive_ratio;
            tmp.energy_total(r)= f.energy_total;

            % 寿命：若未发生死亡则用 T 记为"至少存活到 T"
            tmp.FND(r) = nan_to_T(f.FND, p.T);
            tmp.HND(r) = nan_to_T(f.HND, p.T);
            tmp.LND(r) = nan_to_T(f.LND, p.T);
        end

        for k = 1:numel(metrics)
            x = tmp.(metrics{k});
            res(m).(metrics{k})(vi) = mean(x);
            res(m).([metrics{k} '_std'])(vi) = std(x);
        end
    end
end
end

%% =======================================================================
% 单点重复：用于寿命条形图（FND/HND/LND）
% =======================================================================
function life = sweep_single_point(route_modes, n_runs, base_seed, base)
life = struct();
for m = 1:numel(route_modes)
    mode = route_modes{m};
    FND = zeros(n_runs,1);
    HND = zeros(n_runs,1);
    LND = zeros(n_runs,1);

    for r = 1:n_runs
        p = base;
        p.route_mode = mode;
        p.seed = base_seed + 55500*m + r;
        out = main_baseline(p);
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

% 95% CI
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

hb = bar(ax, Y); %  % 让 MATLAB 自己配色
ng = size(Y,1);
nb = size(Y,2);

% 误差线位置（bar组）
xpos = nan(ng, nb);
for b = 1:nb
    xpos(:,b) = hb(b).XEndPoints;
end
for b = 1:nb
    errorbar(ax, xpos(:,b), Y(:,b), CI(:,b), 'k', 'LineStyle','none', 'CapSize',6);
end

set(ax, 'XTick', 1:ng, 'XTickLabel', modes);
xtickangle(ax, 0);
ylabel(ax, 'Time step');
legend(ax, {'FND','HND','LND'}, 'Location', S.legendLocation, 'Box', S.legendBox, ...
    'FontName', S.fontName, 'FontSize', S.fontSize-1);

tight_inset(ax);
export_figure(gcf, savebase, S);
end

%% =======================================================================
% 工具：NaN 寿命处理
% =======================================================================
function y = nan_to_T(x, T)
if isnan(x), y = T; else, y = x; end
end

%% =======================================================================
% 工具：tight inset 去白边
% =======================================================================
function tight_inset(ax)
outer = ax.OuterPosition;
ti = ax.TightInset;
left = outer(1) + ti(1);
bottom = outer(2) + ti(2);
ax_width = outer(3) - ti(1) - ti(3);
ax_height = outer(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end

%% =======================================================================
% 工具：导出（pdf/eps/png）
% =======================================================================
function export_figure(fig, savebase, S)
set(fig, 'PaperPositionMode','auto');
print(fig, [savebase '.pdf'], '-dpdf', '-vector');
print(fig, [savebase '.eps'], '-depsc', '-vector');
print(fig, [savebase '.png'], ['-r' num2str(S.exportDPI)], '-dpng');
end