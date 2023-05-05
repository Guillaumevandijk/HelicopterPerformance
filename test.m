%test
clear all
[constantParam] = getConstants();
% [v_i, V, Vi_hover] = V_induced(constantParam);
% [P_ideal_act, P_actual_act, P_hov_main] = Hover_power(constantParam);
% [power, V_range, V_endurance] = flight_power(constantParam);


[angles] = Trimmed_condition(constantParam);

m