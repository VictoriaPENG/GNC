function plot_metric_vs_param_abc_paper(resU, resC, resR, x, metric, ylab, xlab, savebase, n_runs, varargin)
% 三拓扑统一 (a)(b)(c) 三联图：uniform/cluster/road

opts = struct('Style',paper_plot_style(),'CI',true,'Legend',true);
opts = parse_opts(opts, varargin{:});
S = opts.Style;

fig = figure('Color','w','Units','centimeters','Position',[2 2 36 12]);
t = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
t.OuterPosition = [0 0.12 1 0.88]; % 底部留给统一 legend

topoTitles = {'(a) Uniform','(b) Cluster','(c) Road'};
RES = {resU, resC, resR};

legend_handles = [];
legend_labels  = {};

for i = 1:3
    ax = nexttile(t,i);
    hold(ax,'on'); grid(ax,'on');
    set(ax,'FontName',S.fontName,'FontSize',S.fontSize,'LineWidth',0.8,'Box','on');
    ax.GridAlpha = S.gridAlpha;

    markers = {'o','s','^','d','v','>','<','p','h','x','+'};

    for m = 1:numel(RES{i})
        y = RES{i}(m).(metric)(:);
        ystd = RES{i}(m).([metric '_std'])(:);

        if opts.CI && n_runs > 1
            ci = S.ciZ * (ystd ./ sqrt(n_runs));
        else
            ci = zeros(size(y));
        end

        mk = markers{1 + mod(m-1, numel(markers))};

        if opts.CI && S.useCIShade && any(ci>0)
            xx = [x(:); flipud(x(:))];
            yy = [y-ci; flipud(y+ci)];
            hfill = fill(ax, xx, yy, 'k', 'FaceAlpha', 0.10, 'EdgeColor','none');
            set(get(get(hfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end

        h = plot(ax, x, y, 'LineWidth', S.lineWidth, 'Marker', mk, ...
            'MarkerSize', S.markerSize, 'DisplayName', RES{i}(m).route_mode);

        if i == 1
            legend_handles(end+1) = h; %#ok<AGROW>
            legend_labels{end+1}  = RES{i}(m).route_mode; %#ok<AGROW>
        end
    end

    title(ax, topoTitles{i}, 'FontWeight','bold');
    xlabel(ax, xlab);
    if i==1, ylabel(ax, ylab); end
end

if opts.Legend
    lgd = legend(t, legend_handles, legend_labels, 'Location','southoutside', 'Orientation','horizontal');
    set(lgd,'Box',S.legendBox,'FontName',S.fontName,'FontSize',S.fontSize-1);
end

paper_plot_style(fig);
export_figure(fig, savebase, S);
end

function opts = parse_opts(opts, varargin)
if mod(numel(varargin),2) ~= 0, error('Optional args must be Name-Value pairs.'); end
for i = 1:2:numel(varargin)
    k = varargin{i}; v = varargin{i+1};
    if isfield(opts,k), opts.(k)=v; else, error('Unknown option: %s', k); end
end
end

function export_figure(fig, savebase, S)
set(fig, 'PaperPositionMode','auto');
print(fig, [savebase '.pdf'], '-dpdf', '-vector');
print(fig, [savebase '.eps'], '-depsc', '-vector');
print(fig, [savebase '.png'], ['-r' num2str(S.exportDPI)], '-dpng');
end
