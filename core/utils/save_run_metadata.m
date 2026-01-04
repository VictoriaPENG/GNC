function save_run_metadata(run_dir, meta, varargin)
%SAVE_RUN_METADATA Save metadata/config for traceable experiments.
%   SAVE_RUN_METADATA(run_dir, meta, 'Config', cfg, 'Results', res)
%
% Saves:
%   - meta.mat (meta + optional cfg/results)
%   - meta.json
%   - README_run.txt (human-readable)

ip = inputParser;
ip.addRequired('run_dir', @(s)ischar(s) || isstring(s));
ip.addRequired('meta', @isstruct);
ip.addParameter('Config', struct(), @isstruct);
ip.addParameter('ResultsSummary', struct(), @isstruct);
ip.addParameter('Extra', struct(), @isstruct);
ip.parse(run_dir, meta, varargin{:});

run_dir = char(ip.Results.run_dir);
cfg = ip.Results.Config;
summary = ip.Results.ResultsSummary;
extra = ip.Results.Extra;

% Save MAT
save(fullfile(run_dir,'meta.mat'), 'meta', 'cfg', 'summary', 'extra');

% Save JSON (best-effort)
json_struct = struct();
json_struct.meta = meta;
json_struct.config = cfg;
json_struct.summary = summary;
json_struct.extra = extra;
try
    txt = jsonencode(json_struct);
    fid = fopen(fullfile(run_dir,'meta.json'), 'w');
    fwrite(fid, txt);
    fclose(fid);
catch
    % ignore JSON errors
end

% Save a human-readable README
fid = fopen(fullfile(run_dir,'README_run.txt'), 'w');
if fid ~= -1
    fprintf(fid, 'Experiment: %s\n', meta.exp_name);
    fprintf(fid, 'Run ID: %s\n', meta.run_id);
    fprintf(fid, 'Timestamp: %s\n', meta.timestamp);
    fprintf(fid, 'Tag: %s\n', meta.tag);
    if isfield(meta,'script') && ~isempty(meta.script)
        fprintf(fid, 'Script: %s\n', meta.script);
    end
    if isfield(meta,'note') && ~isempty(meta.note)
        fprintf(fid, 'Note: %s\n', meta.note);
    end
    fprintf(fid, 'MATLAB: %s (%s)\n', meta.matlab_version, meta.matlab_release);
    fprintf(fid, 'Computer: %s\n', meta.computer);
    if isfield(meta,'hostname') && ~isempty(meta.hostname)
        fprintf(fid, 'Host: %s\n', meta.hostname);
    end
    if isfield(meta,'git') && isfield(meta.git,'available') && meta.git.available
        fprintf(fid, 'Git: %s (branch %s) dirty=%d\n', meta.git.hash, meta.git.branch, meta.git.is_dirty);
    else
        fprintf(fid, 'Git: not available\n');
    end
    fprintf(fid, '\n--- Config ---\n');
    try
        fprintf(fid, '%s\n', evalc('disp(cfg)'));
    catch
        % ignore
    end
    fprintf(fid, '\n--- Summary ---\n');
    try
        fprintf(fid, '%s\n', evalc('disp(summary)'));
    catch
        % ignore
    end
    fclose(fid);
end
end
