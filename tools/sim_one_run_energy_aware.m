function out = sim_one_run_energy_aware(p)
%SIM_ONE_RUN_ENERGY_AWARE Wrapper for the main simulator.
% Separated into an external .m file so PARFOR workers can find it reliably.
out = main_energy_aware_tree(p);
end
