function p = apply_energy_aware_cfg(p)
%APPLY_ENERGY_AWARE_CFG Apply mode-specific configuration for energy-aware tree.

if ~isfield(p,'route_mode') || isempty(p.route_mode)
    return;
end
mode = lower(string(p.route_mode));

if ~isfield(p,'tree_energy_aware') || isempty(p.tree_energy_aware)
    p.tree_energy_aware = false;
end

if mode == "tree_energy"
    p.tree_energy_aware = true;
    if ~isfield(p,'tree_rebuild_period') || isempty(p.tree_rebuild_period) || isinf(p.tree_rebuild_period)
        p.tree_rebuild_period = 50;
    end
    if ~isfield(p,'tree_energy_beta') || isempty(p.tree_energy_beta)
        p.tree_energy_beta = 2.0;
    end
    if ~isfield(p,'tree_energy_min_frac') || isempty(p.tree_energy_min_frac)
        p.tree_energy_min_frac = 0.10;
    end
end
end
