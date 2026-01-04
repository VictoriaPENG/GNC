function root = find_project_root(startPath)
%FIND_PROJECT_ROOT Find project root by searching for startup.m upwards.
%   root = FIND_PROJECT_ROOT(startPath)
%   The project root is defined as the nearest parent folder containing
%   a file named startup.m. If not found, returns startPath.

if nargin < 1 || isempty(startPath)
    startPath = pwd;
end

p = startPath;
root = startPath;
for k = 1:50
    if exist(fullfile(p, 'startup.m'), 'file')
        root = p;
        return;
    end
    parent = fileparts(p);
    if isempty(parent) || strcmp(parent, p)
        break;
    end
    p = parent;
end
end
