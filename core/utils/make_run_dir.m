function [run_dir, meta] = make_run_dir(exp_name, varargin)
%MAKE_RUN_DIR Create a traceable results directory for an experiment.
%   [run_dir, meta] = MAKE_RUN_DIR(exp_name, 'Tag', tag, 'Root', root)
%
% Creates:
%   <root>/results/<exp_name>/<timestamp>__<exp_name>__<tag>
%
% Returns:
%   run_dir: full path to the created directory
%   meta:    struct containing run metadata (timestamp, tag, git hash, etc.)

ip = inputParser;
ip.addRequired('exp_name', @(s)ischar(s) || isstring(s));
ip.addParameter('Tag', 'dev', @(s)ischar(s) || isstring(s));
ip.addParameter('Root', '', @(s)ischar(s) || isstring(s));
ip.addParameter('Script', '', @(s)ischar(s) || isstring(s));
ip.addParameter('Note', '', @(s)ischar(s) || isstring(s));
ip.parse(exp_name, varargin{:});

exp_name = char(ip.Results.exp_name);
tag = char(ip.Results.Tag);
root = char(ip.Results.Root);
script = char(ip.Results.Script);
note = char(ip.Results.Note);

if isempty(root)
    root = find_project_root(pwd);
end

ts = datestr(now, 'yyyymmdd_HHMMSS');
exp_safe = local_sanitize(exp_name);
tag_safe = local_sanitize(tag);
run_id = sprintf('%s__%s__%s', ts, exp_safe, tag_safe);

run_dir = fullfile(root, 'results', exp_safe, run_id);
if ~exist(run_dir, 'dir')
    mkdir(run_dir);
end

meta = struct();
meta.exp_name = exp_name;
meta.tag = tag;
meta.run_id = run_id;
meta.timestamp = ts;
meta.run_dir = run_dir;
meta.root = root;
meta.script = script;
meta.note = note;

% Environment
meta.matlab_version = version;
meta.matlab_release = version('-release');
meta.computer = computer;
meta.hostname = getenv('COMPUTERNAME');
if isempty(meta.hostname)
    meta.hostname = getenv('HOSTNAME');
end

% Git hash (best-effort)
meta.git = local_git_info(root);
end

function s = local_sanitize(s)
if isstring(s), s = char(s); end
s = strrep(s, ' ', '-');
s = regexprep(s, '[^a-zA-Z0-9_\-]', '');
if isempty(s)
    s = 'untitled';
end
end

function g = local_git_info(root)
g = struct('available', false, 'hash', '', 'branch', '', 'is_dirty', false);
try
    [st1, out1] = system(sprintf('cd /d "%s" && git rev-parse HEAD', root));
    if st1 == 0
        g.available = true;
        g.hash = strtrim(out1);
        [st2, out2] = system(sprintf('cd /d "%s" && git rev-parse --abbrev-ref HEAD', root));
        if st2 == 0
            g.branch = strtrim(out2);
        end
        [st3, out3] = system(sprintf('cd /d "%s" && git status --porcelain', root));
        if st3 == 0
            g.is_dirty = ~isempty(strtrim(out3));
        end
    end
catch
    % ignore
end
end
