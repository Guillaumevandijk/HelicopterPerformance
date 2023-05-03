%test
clear all
[Vi_hover, vi_sol, V] = V_induced();
[P_ideal_act, P_actual_act, P_hov_main] = Hover_power();

[constantParam] = getConstants();
tic
[angles] = Trimmed_condition(constantParam);
toc
