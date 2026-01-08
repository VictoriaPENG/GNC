function run_energy_tree()
% ============================================================
% Energy-aware Converge Tree experiment
% ------------------------------------------------------------
% 1) 比较普通汇聚树 (tree) 与能量感知汇聚树 (tree_energy)
% 2) tree_energy 通过 tree_energy_bias 在建树时偏好剩余能量更高的父节点
% 3) 输出结果保存在 results/energy_tree/<timestamp>__energy_tree__<tag> 下
% ============================================================

close all; clc;

% 确保 workers 的 path 包含当前脚本所在目录
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

% ============ 实验配置 ============
tag = getenv('GNC_TAG');
if isempty(tag), tag = 'energy_tree'; end

base = struct();
base.verbose = true;
base.seed = 1;

% 相对 topo_tree 默认配置的差异：缩短仿真步数、适当增加节点数
base.N = 800;
base.T = 1500;
base.route_mode = 'tree';

% 汇聚树参数
base.tree_cost = 'etx';
base.tree_rebuild_period = 250;   % 周期性重建以适配能量耗尽

% 能量感知权重（仅 tree_energy 使用）
energy_bias = 0.6; % 0 表示与 tree 相同，(0,1] 越大越偏向高剩余能量父节点

configs = { ...
    struct('name','tree',         'route_mode','tree',         'tree_energy_bias',0.0), ...
    struct('name','tree_energy',  'route_mode','tree_energy',  'tree_energy_bias',energy_bias) ...
};

% ============ 创建结果目录 ============
[run_dir, meta] = make_run_dir('energy_tree', 'Tag', tag, 'Script', mfilename('fullpath'));
save_run_metadata(run_dir, meta, 'Config', struct('base', base, 'configs', {configs}));

% ============ 运行实验 ============
for k = 1:numel(configs)
    cfg = configs{k};
    fprintf("\\n=== Running %s (route_mode=%s, tree_energy_bias=%.2f) ===\\n", ...
        cfg.name, cfg.route_mode, cfg.tree_energy_bias);

    p = base;
    p.route_mode = cfg.route_mode;
    p.tree_energy_bias = cfg.tree_energy_bias;

    out = main_energy_tree(p);

    save(fullfile(run_dir, sprintf('%s.mat', cfg.name)), 'out', 'p');
end

fprintf("\\nResults saved to: %s\\n", run_dir);
end
