function plot_derived_metric_vs_param_paper(res, x, metricA, metricB, op, ylab, xlab, savebase, n_runs, varargin)
% plot_derived_metric_vs_param_paper
% 派生指标绘图（默认用于 ratio：A/B），并给出 95%CI（误差传播近似）
%
% 例：Energy per delivered ≈ energy_total / PDR
%
% res(m).(metricA), res(m).([metricA '_std'])
% res(m).(metricB), res(m).([metricB '_std'])

opts = struct('Title','', 'YScale','linear','XScale','linear', ...
              'Legend',true,'Style',paper_plot_style(), 'CI',true);
opts = parse_opts(opts, varargin{:});
S = opts.Style;

fig = figure('Color','w','Units','inches', ...
    'Position',[1 1 S.figureW_in S.figureH_in]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on');
set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, ...
    'LineWidth', 0.8, 'Box','on');
if S.useMinorGrid, grid(ax,'minor'); end
ax.GridAlpha = S.gridAlpha; ax.MinorGridAlpha = S.gridAlpha;
set(ax,'XScale',opts.XScale,'YScale',opts.YScale);

markers = {'o','s','^','d','v','>','<','p','h','x','+'};
eps0 = 1e-12;

for m = 1:numel(res)
    A = res(m).(metricA)(:);
    B = res(m).(metricB)(:);
    As = res(m).([metricA '_std'])(:);
    Bs = res(m).([metricB '_std'])(:);

    switch lower(op)
        case 'ratio'
            y = A ./ max(B, eps0);

            % 误差传播（近似）：Var(A/B) ≈ (σA/B)^2 + (AσB/B^2)^2
            ystd = sqrt( (As./max(B,eps0)).^2 + ((A.*Bs)./max(B,eps0).^2).^2 );

        otherwise
            error('Unsupported op: %s (use ''ratio'')', op);
    end

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

    plot(ax, x, y, 'LineWidth', S.lineWidth, 'Marker', mk, ...
        'MarkerSize', S.markerSize, 'DisplayName', res(m).route_mode);

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

if opts.Legend
    lgd = legend(ax, 'Location', S.legendLocation);
    set(lgd, 'Box', S.legendBox, 'FontName', S.fontName, 'FontSize', S.fontSize-1);
end

tight_inset(ax);
paper_plot_style(fig);
export_figure(fig, savebase, S);
end

%% helpers
function opts = parse_opts(opts, varargin)
if mod(numel(varargin),2) ~= 0, error('Optional args must be Name-Value pairs.'); end
for i = 1:2:numel(varargin)
    k = varargin{i}; v = varargin{i+1};
    if isfield(opts,k), opts.(k)=v; else, error('Unknown option: %s', k); end
end
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
