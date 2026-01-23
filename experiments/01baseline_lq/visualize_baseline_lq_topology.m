function [topo, fig] = visualize_baseline_lq_topology(p)
% ============================================================
% visualize_baseline_lq_topology.m
% ------------------------------------------------------------
% 基线实验（baseline_lq）拓扑可视化程序
%
% 使用 generate_topology.m 生成拓扑，并绘制：
%   - 传感器节点（普通/火区）
%   - 基站
%   - 火区范围
%   - cluster/road/hotspot 的辅助信息（可选）
%
% 用法示例：
%   visualize_baseline_lq_topology();  % 默认参数
%   p = struct('topo_mode','cluster','N',500,'show_links',true);
%   visualize_baseline_lq_topology(p);
%   visualize_baseline_lq_topology(struct('topo_mode','all')); % 依次生成三种拓扑图
%
% 输出：
%   topo : generate_topology 生成的结构体
%   fig  : 图句柄
% ============================================================

if nargin < 1
    p = struct();
end
p = fill_defaults(p);

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, '..', '..', 'core'));

if isfield(p, 'topo_modes') && ~isempty(p.topo_modes)
    topo_modes = p.topo_modes;
elseif strcmpi(p.topo_mode, 'all')
    topo_modes = {'uniform','cluster','road'};
else
    topo_modes = {p.topo_mode};
end

topo = cell(numel(topo_modes), 1);
fig = gobjects(numel(topo_modes), 1);

for m = 1:numel(topo_modes)
    p_this = p;
    p_this.topo_mode = topo_modes{m};
    if ~isempty(p.save_path) && numel(topo_modes) > 1
        [save_dir, save_name, save_ext] = fileparts(p.save_path);
        if isempty(save_ext)
            save_ext = '.png';
        end
        p_this.save_path = fullfile(save_dir, sprintf('%s_%s%s', save_name, p_this.topo_mode, save_ext));
    end

    [topo{m}, fig(m)] = plot_single_topology(p_this);
end

if isscalar(topo_modes)
    topo = topo{1};
end

end

function h = draw_circle(center, radius, style, width)
theta = linspace(0, 2*pi, 200);
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);
h = plot(x, y, style, 'LineWidth', width);
end

function [topo, fig] = plot_single_topology(p)
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

pos = topo.pos_sensors;
bs_pos = topo.BS_pos;

is_fire_node = vecnorm(pos - p.fire_pos, 2, 2) <= p.Rf;

fig = figure(...
    'Units','centimeters', ...
    'Position',[2 2 20 20], ...
    'Color','w');
hold on; axis equal;
axis([0 p.Lx 0 p.Ly]);
xlim([0 p.Lx]); ylim([0 p.Ly]);
title(sprintf('baseline\\_lq Topology (%s)', p.topo_mode), 'Interpreter','none');
xlabel('X (m)'); ylabel('Y (m)');

legend_handles = [];
legend_labels = {};

if p.show_links
    pos_all = topo.pos_all;
    D = squareform(pdist(pos_all));
    Adj = (D <= p.Rc) & (D > 0);
    [i_idx, j_idx] = find(triu(Adj(1:p.N, 1:p.N)));
    for k = 1:numel(i_idx)
        i = i_idx(k); j = j_idx(k);
        plot([pos(i,1) pos(j,1)], [pos(i,2) pos(j,2)], ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    end
    legend_handles(end+1) = plot(nan, nan, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5); 
    legend_labels{end+1} = 'Links'; 


end

h_nodes = scatter(pos(~is_fire_node,1), pos(~is_fire_node,2), ...
    p.node_size, p.node_color, 'filled', 'MarkerEdgeColor','none');
h_fire_nodes = scatter(pos(is_fire_node,1), pos(is_fire_node,2), ...
    p.node_size + 4, p.fire_node_color, 'filled', 'MarkerEdgeColor','none');

% ✅ 半透明（R2025b 支持 MarkerFaceAlpha）
if isfield(p,'node_alpha') && ~isempty(p.node_alpha)
    h_nodes.MarkerFaceAlpha = p.node_alpha;
    h_fire_nodes.MarkerFaceAlpha = min(1, p.node_alpha + 0.2);
end


h_bs = plot(bs_pos(1), bs_pos(2), 'p', 'MarkerSize', 14, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', p.bs_color);

h_fire = draw_circle(p.fire_pos, p.Rf, p.fire_circle_style, p.fire_circle_width);

if isfield(topo.meta,'cluster_centers')
    centers = topo.meta.cluster_centers;
    plot(centers(:,1), centers(:,2), 'kx', 'LineWidth', 1.5, 'MarkerSize', 8);
end

if isfield(topo.meta,'roads')
    roads = topo.meta.roads;
    for r = 1:numel(roads)
        poly = roads{r};
        plot(poly(:,1), poly(:,2), '-', 'Color', p.road_color, 'LineWidth', 1.5);
    end
end

if isfield(topo.meta,'hotspot_center')
    hc = topo.meta.hotspot_center;
    plot(hc(1), hc(2), 'o', 'MarkerSize', 10, ...
        'MarkerEdgeColor', p.hotspot_color, 'LineWidth', 1.5);
end

legend_handles = [legend_handles, h_nodes, h_fire_nodes, h_bs, h_fire];
legend_labels = [legend_labels, {'Sensors','Fire sensors','Base station','Fire region'}];
lgd = legend(legend_handles, legend_labels, ...
    'Location','southoutside', 'Orientation','horizontal');
lgd.Box = 'off';

grid on;

if ~isempty(p.save_path)
    [save_dir,~,~] = fileparts(p.save_path);
    if ~isempty(save_dir) && ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    exportgraphics(fig, p.save_path, 'Resolution', 300);  % 300dpi 论文常用
end
end

function p = fill_defaults(p)
def = default_system_params();

fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(p, f)
        p.(f) = def.(f);
    end
end
end
