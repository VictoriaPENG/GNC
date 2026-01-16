function plot_model_validation_figs(p)
% ============================================================
% plot_model_validation_figs.m
% ------------------------------------------------------------
% 基线实验（baseline_lq）建模验证专用画图程序
% 输出图：
%   F1-1 拓扑可视化
%   F1-2 连通率 vs Rc
%   F1-3 PRR/ETX vs 距离
%   F1-4 队列溢出 vs 负载
%   F1-5 能耗分解 / 累计能耗
%
% 保存图片到本地，默认白色背景。
% ============================================================


if nargin < 1
    p = struct();
end
p = fill_defaults(p);

thisDir = fileparts(mfilename('fullpath'));
%% Path setup (robust & project-structure aware)
this_dir  = fileparts(mfilename('fullpath'));  % ...\GNC\plots
proj_root = fileparts(this_dir);               % ...\GNC
p.outdir  = fullfile(proj_root, 'results', 'figs_model_validation');
if ~exist(p.outdir,'dir'), mkdir(p.outdir); end
fprintf('[plot_model_validation_figs] outdir = %s\n', p.outdir);



%% F1-1 拓扑可视化
fig1_path = fullfile(p.outdir, 'F1-1_topology.png');
apply_paper_plot_style;

visualize_baseline_lq_topology(struct( ...
    'topo_mode', p.topo_mode, ...
    'N', p.N, ...
    'Lx', p.Lx, ...
    'Ly', p.Ly, ...
    'Rc', p.Rc, ...
    'BS_pos', p.BS_pos, ...
    'fire_pos', p.fire_pos, ...
    'Rf', p.Rf, ...
    'seed', p.seed, ...
    'show_links', p.show_links, ...
    'save_path', fig1_path));

%% F1-2 连通率 vs Rc
Rc_list = p.Rc_list;
conn_ratio = zeros(numel(Rc_list),1);
for k = 1:numel(Rc_list)
    conn_ratio(k) = connectivity_ratio(p, Rc_list(k));
end
fig = figure('Color','w');
h = plot(Rc_list, conn_ratio, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
xlabel('R_c (m)');
ylabel('Connectivity ratio');
title('F1-2 Connectivity vs R_c');
grid on;
apply_paper_plot_style;

save_figure(fig, fullfile(p.outdir, 'F1-2_connectivity_vs_Rc'));

%% F1-3 PRR / ETX vs 距离
d = linspace(0, p.prr_dist_max, 200);
prr = exp(-p.alpha * d);
prr = max(prr, p.prr_floor);
prr = min(prr, 1);
etx = 1 ./ max(prr, 1e-6);

fig = figure('Color','w');
yyaxis left;
h1 = plot(d, prr, '-', 'LineWidth', 1.8);
ylabel('PRR');
ylim([0 1.05]);
yyaxis right;
h2 = plot(d, etx, '--', 'LineWidth', 1.8);
ylabel('ETX');
xlabel('Distance (m)');
title('F1-3 PRR / ETX vs Distance');
grid on;
legend([h1 h2], {'PRR','ETX'}, 'Location','best');
apply_paper_plot_style;

save_figure(fig, fullfile(p.outdir, 'F1-3_prr_etx_vs_distance'));

%% F1-4 队列溢出 vs 负载
lambda_list = p.lambda_list;
q_overflow_mean = zeros(numel(lambda_list),1);
for k = 1:numel(lambda_list)
    overflows = zeros(p.n_runs,1);
    for r = 1:p.n_runs
        p_run = p;
        p_run.lambda = lambda_list(k);
        p_run.seed = p.seed + r - 1;
        p_run.energy_trace = false;
        out = main_baseline(p_run);
        overflows(r) = out.diag.q_overflow;
    end
    q_overflow_mean(k) = mean(overflows);
end
fig = figure('Color','w');
plot(lambda_list, q_overflow_mean, '-s', 'LineWidth', 1.8, 'MarkerSize', 6);
xlabel('\lambda (packet generation probability)');
ylabel('Queue overflow count');
title('F1-4 Queue Overflow vs Load');
grid on;
apply_paper_plot_style;

save_figure(fig, fullfile(p.outdir, 'F1-4_queue_overflow_vs_load'));

%% F1-5 能耗分解 / 累计能耗
p_energy = p;
p_energy.energy_trace = true;
out = main_baseline(p_energy);

fig = figure('Color','w');
t = (1:p_energy.T)';

% 左图：能耗分解（堆叠柱状图）
subplot(1,2,1);
idle = out.final.energy_idle;
tx   = out.final.energy_tx;
rx   = out.final.energy_rx;

Y = [idle, tx, rx];   % 1×3
x = 1;                % 只有一组：Total
hb = bar(x, Y, 'stacked');  % 返回 1×3 Bar 对象，legend 不会再多余/缺失

xticks(1);
xticklabels({'Total'});
xlim([0.5 1.5]);
legend(hb(:), {'Idle','Tx','Rx'}, 'Location','best');

ylabel('Energy');
title('F1-5 Energy Breakdown');
grid on;

% 右图：累计能耗随时间
subplot(1,2,2);
plot(t, out.time.energy.cum, 'LineWidth', 1.8);
xlabel('Time step');
ylabel('Cumulative Energy');
title('F1-5 Cumulative Energy');
grid on;

apply_paper_plot_style;

save_figure(fig, fullfile(p.outdir, 'F1-5_energy_breakdown_cumulative'));


end

function ratio = connectivity_ratio(p, Rc)
cfgTopo = struct();
cfgTopo.N = p.N;
cfgTopo.Lx = p.Lx;
cfgTopo.Ly = p.Ly;
cfgTopo.seed = p.seed;
cfgTopo.topo_mode = p.topo_mode;
cfgTopo.BS_pos = p.BS_pos;
cfgTopo.fire_pos = p.fire_pos;
cfgTopo.Rf = p.Rf;

cfgTopo.hotspot_enable = p.hotspot_enable;
cfgTopo.hotspot_ratio  = p.hotspot_ratio;
cfgTopo.hotspot_sigma  = p.hotspot_sigma;
cfgTopo.hotspot_center = p.hotspot_center;

cfgTopo.Kc = p.Kc;
cfgTopo.cluster_sigma = p.cluster_sigma;
cfgTopo.cluster_center_mode = p.cluster_center_mode;
cfgTopo.cluster_mix_uniform_ratio = p.cluster_mix_uniform_ratio;

cfgTopo.road_count = p.road_count;
cfgTopo.road_width = p.road_width;
cfgTopo.road_mix_uniform_ratio = p.road_mix_uniform_ratio;
cfgTopo.road_node_ratio = p.road_node_ratio;
if isfield(p,'roads') && ~isempty(p.roads)
    cfgTopo.roads = p.roads;
end

topo = generate_topology(cfgTopo);
pos_all = topo.pos_all;
bs_idx = p.N + 1;
dist = squareform(pdist(pos_all));
Adj = (dist <= Rc) & (dist > 0);
G = graph(Adj);
dist_hop = distances(G, 1:p.N, bs_idx);
ratio = mean(isfinite(dist_hop));
end

function apply_paper_style(fig)
% 统一设置：论文级字体/线宽/背景，并移除坐标区工具栏（避免 export 时把 UI 导出到图片里）
if nargin < 1 || isempty(fig)
    fig = gcf;
end

set(fig, 'Color', 'w');

% 统一全局排版
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    ax = axs(i);
    % 白底 + 边框 + 网格更清爽
    ax.Box = 'on';
    ax.LineWidth = 1.0;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 11;
    ax.LabelFontSizeMultiplier = 1.1;
    ax.TitleFontSizeMultiplier = 1.1;

    % 移除 Axes Toolbar（R2020b+）
    try
        ax.Toolbar = [];
    catch
        % older MATLAB: ignore
    end
    try
        ax.Toolbar.Visible = 'off';
    catch
        % some versions do not expose Visible
    end

    % 网格线更适合论文（细一点）
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.10;
end

% 统一线条/标记风格（已画好的对象也一起改）
lns = findall(fig, 'Type', 'line');
for i = 1:numel(lns)
    if lns(i).LineWidth < 1.5
        lns(i).LineWidth = 1.8;
    end
    if ~strcmp(lns(i).Marker, 'none') && lns(i).MarkerSize < 6
        lns(i).MarkerSize = 6;
    end
end

% 图例排版
legs = findall(fig, 'Type', 'legend');
for i = 1:numel(legs)
    legs(i).Box = 'off';
    legs(i).FontName = 'Times New Roman';
    legs(i).FontSize = 10;
end

% 确保向量导出质量（pdf/eps）
try
    set(fig, 'Renderer', 'painters');
catch
end
end

function save_figure(fig, path)
% 论文级导出：
% 1) 默认 300dpi png
% 2) 同名导出一份 pdf（矢量）
% 3) 自动移除坐标区工具栏，避免导出 UI
apply_paper_style(fig);

[save_dir, base, ext] = fileparts(path);
if ~isempty(save_dir) && ~exist(save_dir,'dir')
    mkdir(save_dir);
end

if isempty(ext)
    ext = '.png';
end

% png（位图，适合论文插图）
png_path = fullfile(save_dir, [base '.png']);
try
    exportgraphics(fig, png_path, 'Resolution', 300, 'BackgroundColor', 'white');
catch
    % 兼容老版本
    print(fig, png_path, '-dpng', '-r300');
end

% pdf（矢量，适合 LaTeX / Word）
pdf_path = fullfile(save_dir, [base '.pdf']);
try
    exportgraphics(fig, pdf_path, 'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    % 兼容老版本
    print(fig, pdf_path, '-dpdf', '-vector');
end

close(fig);
end

function p = fill_defaults(p)
def = struct();
def.Lx = 1000; def.Ly = 1000;
def.N  = 700;
def.Rc = 70;
def.T  = 2000;

def.BS_pos   = [500, 500];
def.fire_pos = [250, 250];
def.Rf = 220;

def.lambda = 0.05;
def.alpha = 0.01;
def.prr_floor = 0.08;

def.seed = 1;
def.route_mode = 'mix';
def.verbose = false;

def.topo_mode = 'uniform';

def.hotspot_enable = false;
def.hotspot_ratio  = 0.25;
def.hotspot_sigma  = min(def.Lx,def.Ly)*0.05;
def.hotspot_center = def.fire_pos;

def.Kc = 6;
def.cluster_sigma = min(def.Lx,def.Ly)*0.06;
def.cluster_center_mode = 'random';
def.cluster_mix_uniform_ratio = 0.2;

def.road_count = 2;
def.road_width = min(def.Lx,def.Ly)*0.01;
def.road_mix_uniform_ratio = 0.2;
def.road_node_ratio = 0.8;

def.show_links = false;
def.outdir = fullfile(pwd, 'figs_model_validation');
def.Rc_list =  [20 30 40 50 60 80 100 120 150 180 220];
def.lambda_list = [0.01 0.02 0.05 0.08 0.10];
def.prr_dist_max = 300;
def.n_runs = 5;

fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(p, f)
        p.(f) = def.(f);
    end
end
end
