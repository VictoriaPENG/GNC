function plot_metric_vs_param_paper(res, x, metric, ylab, xlab, savebase, n_runs, varargin)
% 论文级对比曲线绘图（均值 + 95%CI 阴影带/误差棒）并导出 PDF/EPS/PNG
%
% res(m).(metric)      : mean vector
% res(m).([metric '_std']) : std vector
% n_runs               : 每个点重复次数，用于 CI
%
% 可选参数（Name-Value）：
% 'Title'      : 图标题（默认不写，论文通常靠 caption）
% 'YScale'     : 'linear' or 'log'
% 'XScale'     : 'linear' or 'log'
% 'Legend'     : true/false
% 'Style'      : style struct（默认 paper_plot_style()）
% 'CI'         : true/false（是否显示CI）
%
% savebase 不带后缀，例如 'figs/PDR_vs_lambda'

opts = struct('Title','', 'YScale','linear','XScale','linear', ...
              'Legend',true,'Style',paper_plot_style(), 'CI',true);
opts = parse_opts(opts, varargin{:});
S = opts.Style;

% Figure
fig = figure('Color','w','Units','inches', ...
    'Position',[1 1 S.figureW_in S.figureH_in]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on');

set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, ...
    'LineWidth', 0.8, 'Box','on');

if S.useMinorGrid
    grid(ax,'minor');
end
ax.GridAlpha = S.gridAlpha;
ax.MinorGridAlpha = S.gridAlpha;

set(ax,'XScale',opts.XScale,'YScale',opts.YScale);

% 为可复现，固定 marker 序列（不指定具体颜色，让 MATLAB 自动循环）
markers = {'o','s','^','d','v','>','<','p','h','x','+'};

for m = 1:numel(res)
    y = res(m).(metric)(:);
    ystd = res(m).([metric '_std'])(:);

    % 95% CI
    if opts.CI && n_runs > 1
        ci = S.ciZ * (ystd ./ sqrt(n_runs));
    else
        ci = zeros(size(y));
    end

    mk = markers{1 + mod(m-1, numel(markers))};

    % 先画 CI 阴影带（更论文感）
    if opts.CI && S.useCIShade && any(ci>0)
        xx = [x(:); flipud(x(:))];
        yy = [y-ci; flipud(y+ci)];
        hfill = fill(ax, xx, yy, 'k', 'FaceAlpha', 0.10, 'EdgeColor','none');
        % 不让阴影进 legend
        set(get(get(hfill,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end

    % 画均值曲线
    plot(ax, x, y, 'LineWidth', S.lineWidth, 'Marker', mk, ...
        'MarkerSize', S.markerSize, 'DisplayName', res(m).route_mode);

    % 可选：误差棒（如果你更喜欢 errorbar 风格）
    if opts.CI && ~S.useCIShade && any(ci>0)
        er = errorbar(ax, x, y, ci, 'LineStyle','none', 'CapSize',6);
        set(get(get(er,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

xlabel(ax, xlab, 'FontName', S.fontName, 'FontSize', S.fontSize);
ylabel(ax, ylab, 'FontName', S.fontName, 'FontSize', S.fontSize);

if ~isempty(opts.Title)
    title(ax, opts.Title, 'FontName', S.fontName, 'FontSize', S.titleSize, 'FontWeight','normal');
end

% Legend
if opts.Legend
    lgd = legend(ax, 'Location', S.legendLocation);
    set(lgd, 'Box', S.legendBox, 'FontName', S.fontName, 'FontSize', S.fontSize-1);
end

% 去白边（论文排版更干净）
tight_inset(ax);

% 导出（矢量优先）
export_figure(fig, savebase, S);

end

%% -------- helpers --------
function opts = parse_opts(opts, varargin)
if mod(numel(varargin),2) ~= 0
    error('Optional args must be Name-Value pairs.');
end
for i = 1:2:numel(varargin)
    k = varargin{i};
    v = varargin{i+1};
    if isfield(opts, k)
        opts.(k) = v;
    else
        error('Unknown option: %s', k);
    end
end
end

function tight_inset(ax)
% 让 axes 更紧凑，减少留白（适合论文排版）
outer = ax.OuterPosition;
ti = ax.TightInset;
left = outer(1) + ti(1);
bottom = outer(2) + ti(2);
ax_width = outer(3) - ti(1) - ti(3);
ax_height = outer(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end

function export_figure(fig, savebase, S)
% 矢量 PDF（首选）、EPS（期刊常要）、PNG（预览）
set(fig, 'PaperPositionMode','auto');

% PDF (vector)
print(fig, [savebase '.pdf'], '-dpdf', '-painters');

% EPS (vector)
print(fig, [savebase '.eps'], '-depsc', '-painters');

% PNG (raster preview)
print(fig, [savebase '.png'], ['-r' num2str(S.exportDPI)], '-dpng');
end
