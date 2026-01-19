function S = paper_plot_style(fig)
%PAPER_PLOT_STYLE Paper unified style + (optional) apply to figure.
%
% Usage:
%   S = paper_plot_style();        % only get style struct
%   S = paper_plot_style(gcf);     % get style struct AND apply to current figure

% =========================
% 1) Style config (your S)
% =========================
S.fontName   = 'Times New Roman';
S.fontSize   = 11;
S.titleSize  = 11;

S.lineWidth  = 1.6;
S.markerSize = 6;
S.gridAlpha  = 0.15;

S.figureW_in = 3.5;     % single column, inches (double column ~7.2)
S.figureH_in = 2.6;

S.legendLocation = 'best';
S.legendBox      = 'off';

S.useMinorGrid = true;
S.useCIShade   = true;
S.ciZ          = 1.96;

S.exportDPI    = 300;

% =========================
% 2) Apply (merged from apply_paper_plot_style)
% =========================
if nargin < 1 || isempty(fig) || ~ishandle(fig)
    return; % no figure handle => only return S
end

set(fig, 'Color','w', 'Units','inches', ...
    'Position',[1 1 S.figureW_in S.figureH_in]);

axs = findall(fig, 'Type','axes');
for k = 1:numel(axs)
    ax = axs(k);

    set(ax, 'FontName', S.fontName, 'FontSize', S.fontSize, ...
        'LineWidth', 1.0, 'Box','on', 'Color','w');

    grid(ax,'on');
    ax.GridAlpha = S.gridAlpha;

    if isfield(S,'useMinorGrid') && S.useMinorGrid
        grid(ax,'minor');
    end

    % title/labels
    t  = get(ax,'Title');  if ~isempty(t),  set(t, 'FontName',S.fontName, 'FontSize',S.titleSize); end
    xl = get(ax,'XLabel'); if ~isempty(xl), set(xl,'FontName',S.fontName, 'FontSize',S.fontSize); end
    yl = get(ax,'YLabel'); if ~isempty(yl), set(yl,'FontName',S.fontName, 'FontSize',S.fontSize); end
end

ln = findall(fig, 'Type','line');
for k = 1:numel(ln)
    if isprop(ln(k),'LineWidth'),  ln(k).LineWidth  = S.lineWidth;  end
    if isprop(ln(k),'MarkerSize'), ln(k).MarkerSize = S.markerSize; end
end

lgd = findall(fig, 'Type','legend');
for k = 1:numel(lgd)
    set(lgd(k), 'Location', S.legendLocation, 'Box', S.legendBox, ...
        'FontName', S.fontName, 'FontSize', S.fontSize);
end

drawnow;
end
