function apply_paper_plot_style()
% 全局论文级绘图风格（startup 调用）

set(groot, 'defaultFigureColor', 'w');

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1);

set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 8);

set(groot, 'defaultAxesGridLineStyle', '--');
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');

set(groot, 'defaultLegendBox', 'off');
set(groot, 'defaultLegendFontSize', 11);
set(groot, 'defaultLegendTextColor', 'k');

disp('[startup] Paper plot style applied.');
end

