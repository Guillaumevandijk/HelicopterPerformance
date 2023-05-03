%test
clear all
[constantParam] = getConstants();
%[v_i, V] = V_induced(constantParam);
%[P_ideal_act, P_actual_act, P_hov_main] = Hover_power(constantParam);
[power] = flight_power(constantParam);

%tic
%[angles] = Trimmed_condition(constantParam);
%toc
