function fig = plot_topology_abc_paper(p)
% ============================================================
% Figure 1(a)(b)(c) 模板：uniform / cluster / road
% - 不影响其他图：只在这里控制 figure 大小、marker、布局
% - 统一坐标范围、比例、字体、子图标注
% ============================================================

if nargin < 1, p = struct(); end
p = fill_defaults_local(p);

modes = {'uniform','cluster','road'};
mode_titles = {'Uniform','Cluster','Road'};

% --- 大画布：三联图（横向），拓扑图必须更大 ---
fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 36 13]);
t = tiledlayout(fig, 1, 3, 'Padding','compact', 'TileSpacing','compact');
t.OuterPosition = [0 0.10 1 0.90];   % ✅ 底部预留 10% 放 legend

% 为了统一图例：我们只在最后统一 legend（放到底部/外侧）
legend_handles = [];
legend_labels  = {};

for i = 1:3
    ax = nexttile(t, i);
    hold(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [0 p.Lx]); ylim(ax, [0 p.Ly]);
    grid(ax, 'on');
    set(ax, 'FontSize', 14, 'LineWidth', 1.0);
    title(ax, mode_titles{i}, 'FontSize', 18, 'FontWeight','bold');


    % --- 生成对应拓扑 ---
    cfgTopo = struct();
    cfgTopo.N = p.N;
    cfgTopo.Lx = p.Lx;
    cfgTopo.Ly = p.Ly;
    cfgTopo.seed = p.seed + i;     % ✅ 三幅图避免完全一样（也可去掉）
    cfgTopo.topo_mode = modes{i};
    cfgTopo.BS_pos = p.BS_pos;
    cfgTopo.fire_pos = p.fire_pos;
    cfgTopo.Rf = p.Rf;

    % cluster / road / hotspot 参数（保持与你现有一致）
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

    % --- 画 links（可选，默认关，避免太乱）---
    if p.show_links
        pos_all = topo.pos_all;
        D = squareform(pdist(pos_all));
        Adj = (D <= p.Rc) & (D > 0);
        [ii, jj] = find(triu(Adj(1:p.N, 1:p.N)));
        for k = 1:numel(ii)
            plot(ax, [pos(ii(k),1) pos(jj(k),1)], [pos(ii(k),2) pos(jj(k),2)], ...
                'Color', [0.85 0.85 0.85], 'LineWidth', 0.3);
        end
    end

    % --- 点太密：拓扑专用 marker（比指标图更小）---
    node_size = p.node_size_abc;  % ✅ 单独给三联图用的小点
    h_nodes = scatter(ax, pos(~is_fire_node,1), pos(~is_fire_node,2), ...
        node_size, p.node_color, 'filled');
    h_fire_nodes = scatter(ax, pos(is_fire_node,1), pos(is_fire_node,2), ...
        node_size+6, p.fire_node_color, 'filled');

    h_bs = plot(ax, bs_pos(1), bs_pos(2), 'p', 'MarkerSize', 11, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', p.bs_color);

    % fire region circle
    draw_circle(ax, p.fire_pos, p.Rf, p.fire_circle_style, p.fire_circle_width);

    % cluster centers / roads / hotspot 等辅助信息
    if isfield(topo.meta,'cluster_centers')
        centers = topo.meta.cluster_centers;
        plot(ax, centers(:,1), centers(:,2), 'kx', 'LineWidth', 1.2, 'MarkerSize', 6);
    end
    if isfield(topo.meta,'roads')
        roads = topo.meta.roads;
        for r = 1:numel(roads)
            poly = roads{r};
            plot(ax, poly(:,1), poly(:,2), '-', 'Color', p.road_color, 'LineWidth', 1.2);
        end
    end
    if isfield(topo.meta,'hotspot_center')
        hc = topo.meta.hotspot_center;
        plot(ax, hc(1), hc(2), 'o', 'MarkerSize', 8, ...
            'MarkerEdgeColor', p.hotspot_color, 'LineWidth', 1.2);
    end

    % --- 标题 & 子图标注 (a)(b)(c) ---
    title(ax, mode_titles{i}, 'Interpreter','none');
    text(ax, 0.02*p.Lx, 0.96*p.Ly, sprintf('(%c)', 'a'+(i-1)), ...
        'FontWeight','bold', 'FontSize', 14, 'VerticalAlignment','top');

    % 统一坐标标签：只在左下角显示，避免挤
    if i == 1
        ylabel(ax, 'Y (m)');
    else
        ax.YTickLabel = [];
    end
    xlabel(ax, 'X (m)');

    % --- 收集图例句柄（只收一次）---
    if i == 1
        legend_handles = [h_nodes, h_fire_nodes, h_bs];
        legend_labels  = {'Sensors','Fire sensors','Base station'};
    end
end

% 统一图例：放在底部横排（不遮挡）
ax1 = findobj(fig, 'Type', 'axes'); 
ax1 = ax1(end);

lgd = legend(ax1, legend_handles, legend_labels, ...
    'Orientation','horizontal');
lgd.Box = 'off';

% ✅ 把 legend 放进 tiledlayout 的 south tile（不挤子图）
try
    lgd.Layout.Tile = 'south';
    lgd.Layout.TileSpan = [1 3];   % 横跨 3 列
catch
    % 老版本不支持 Layout.Tile，就退回到 southoutside
    lgd.Location = 'southoutside';
end


% 保存（建议 exportgraphics，尺寸更稳定）
if ~isempty(p.save_path)
    [save_dir,~,~] = fileparts(p.save_path);
    if ~isempty(save_dir) && ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    exportgraphics(fig, p.save_path, 'Resolution', 400);
end

end

% ======== helpers ========
function draw_circle(ax, center, radius, style, width)
theta = linspace(0, 2*pi, 200);
x = center(1) + radius*cos(theta);
y = center(2) + radius*sin(theta);
plot(ax, x, y, style, 'LineWidth', width);
end

function p = fill_defaults_local(p)
def = default_system_params();

% 拓扑参数（沿用你工程默认）
def.topo_mode = 'uniform';
def.show_links = false;

% hotspot
def.hotspot_enable = false;
def.hotspot_ratio  = 0.25;
def.hotspot_sigma  = min(def.Lx,def.Ly)*0.05;
def.hotspot_center = def.fire_pos;

% cluster
def.Kc = 6;
def.cluster_sigma = min(def.Lx,def.Ly)*0.06;
def.cluster_center_mode = 'random';
def.cluster_mix_uniform_ratio = 0.2;

% road
def.road_count = 2;
def.road_width = min(def.Lx,def.Ly)*0.01;
def.road_mix_uniform_ratio = 0.2;
def.road_node_ratio = 0.8;

% 颜色/视觉
def.node_color = [0.2 0.2 0.2];
def.fire_node_color = [0.85 0.2 0.2];
def.bs_color = [0.1 0.4 0.9];
def.road_color = [0.2 0.6 0.2];
def.hotspot_color = [0.7 0.2 0.7];
def.fire_circle_style = '--r';
def.fire_circle_width = 1.3;

% ✅ 关键：三联图专用 marker（不影响你别的图）
def.node_size_abc = 10;

% 输出
def.save_path = '';

fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(p, fn{k})
        p.(fn{k}) = def.(fn{k});
    end
end
end
