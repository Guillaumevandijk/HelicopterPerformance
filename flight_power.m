function [power] = flight_power(constantParam)
%This function calculates the total required power in forward flight using
%BEM
[v_in, V] = V_induced(constantParam);

%Rotor parasite power and profile drag power
blade_sol = constantParam.sigma;
C_Dp = constantParam.C_dp_main;
rho = constantParam.rho;
omega_main = constantParam.Omega;
R_main = constantParam.R_tip;
A_fus = constantParam.S_eq;
W = constantParam.mass*9.81;
Tail_arm = constantParam.Tailarm;

power = zeros(1, (length(V)));

for count = 1:length(V)


%Induced power
k = constantParam.k_main;
Drag_fus = 0.5 * rho * V(count)^2 * A_fus; %Uses equivalent flat plate area
T = sqrt(W^2 + Drag_fus^2);
v_i = v_in(count);
P_i = k*T*v_i;


mu = V(count)/(omega_main*R_main);

P_par = blade_sol*(C_Dp/8)*rho*((omega_main*R_main)^3)*pi*(R_main^2)*(1+4.65*mu^2);

%Fuselage drag power Not needed according to assignment?
P_d = Drag_fus*V(count);

P_tot_main = P_i + P_par + P_d;

%----tail rotor-----------------
M_tail = P_tot_main/omega_main;
T_tail = M_tail/Tail_arm;
R_tail = constantParam.R_tail;
vi_tail = sqrt(T_tail/(2*rho*pi*(R_tail^2)));
k_tail = constantParam.k_tail;

%Induced tail power
P_i_tail = k_tail*T_tail*vi_tail;

Solidity_tail = constantParam.sigma_tail;
C_Dp_tail = constantParam.C_dp_tail;
Omega_tail = constantParam.Omega_tail;
mu_tail = 0; %Dit is 0 toch? Geen yawing motion of sideslip dus advance ratio = 0

P_par_tail = Solidity_tail*(C_Dp_tail/8)*rho*((Omega_tail*R_tail)^3)*pi*(R_tail^2)*(1+4.65*mu_tail^2);

P_tot_tail = P_i_tail + P_par_tail;

P_tot_all = P_tot_main + P_tot_tail;

power(count) = P_tot_all;

end

figure
%plot(V, v_i, 'DisplayName', 'Induced velocity', LineWidth=1.0)
plot(V, power, LineWidth=1.0)
hold on
%legend('Interpreter','latex')
grid
title('Required total power in forward flight', 'Interpreter','latex')
xlabel('Forward velocity [m/s]')
%xlim([0 100])
ylabel('Required power [W]')
end