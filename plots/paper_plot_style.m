function S = paper_plot_style()
% 论文规范统一风格（IEEE/Sensors 常用）
S.fontName   = 'Times New Roman';
S.fontSize   = 11;     % 正文图常用 9~11
S.titleSize  = 11;
S.lineWidth  = 1.6;
S.markerSize = 6;
S.gridAlpha  = 0.15;

S.figureW_in = 3.5;    % 单栏图宽度（inches），双栏可用 7.2
S.figureH_in = 2.6;

S.legendLocation = 'best';
S.legendBox = 'on';

S.useMinorGrid = true;
S.useCIShade   = true;    % true=阴影置信带；false=误差棒
S.ciZ          = 1.96;    % 95% normal approx

% 导出
S.exportDPI = 300;

end
