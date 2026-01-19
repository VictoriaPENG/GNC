function save_figure(fig, filepath_no_ext)
%SAVE_FIGURE Paper-quality export: PNG(300dpi)+PDF(vector)
% Usage:
%   save_figure(gcf, fullfile(outdir, 'F1-5_energy_breakdown'));

    if nargin < 2 || isempty(filepath_no_ext)
        error('save_figure: filepath_no_ext is required');
    end

    outdir = fileparts(filepath_no_ext);
    if ~isempty(outdir) && ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % White background for paper
    set(fig, 'Color', 'w');

    exportgraphics(fig, [filepath_no_ext '.png'], 'Resolution', 300);
    exportgraphics(fig, [filepath_no_ext '.pdf'], 'ContentType', 'vector');
end
