function [P_ideal_act, P_actual_act, P_hov_main] = Hover_power(constantParam)
%This function calculates the ideal and actual power for hover

W = constantParam.mass*9.81;
rho = constantParam.rho;
R = constantParam.R_tip;
FM = constantParam.FM;

%------------------ACT-------------------%
%Hover power with ACT
P_ideal_act = W*sqrt(W/(2*pi*rho*R^2));
P_actual_act = P_ideal_act/FM;
%------------------ACT-------------------%

%------------------BEM-------------------%
%Hover power with BEM
k_main = constantParam.k_main;
T_main = W; % T equals weight
vi_main = sqrt(W/(2*rho*pi*(R^2)));
R_main = R;
Omega_main = constantParam.Omega; %rps
C_dp_main = 0.0057;
blade_sol_main = constantParam.sigma;

P_i_main = k_main*T_main*vi_main;
P_p_main = blade_sol_main*C_dp_main/8*rho*((Omega_main*R_main)^3)*pi*R_main^2;
P_hov_main = P_i_main+P_p_main;

%------------------BEM-------------------%

%P_i_tail = k_tail*T_tail*vi_tail;
%P_p_tail = blade_sol_tail*C_dp_tail/8*rho*((Omega_tail*R_tail)^3)*pi*R_tail^2;
%P_hov_tail = p_i_mtail+P_p_tail;
end