function startup()
% GNC 工程统一启动：初始化路径 + 清缓存
root = fileparts(mfilename('fullpath'));
cd(root);

addpath(genpath(fullfile(root,'core')));
addpath(genpath(fullfile(root,'plots')));
addpath(genpath(fullfile(root,'experiments')));
addpath(genpath(fullfile(root,'tools')));

rehash;
disp("GNC project initialized: core/ plots/ experiments/ added to path");
end
