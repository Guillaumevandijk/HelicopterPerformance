function [v_i, V] = V_induced(constantParam)
%This function calculates the induced velocity in hover and forward flight
%using the Actuator Disc Theory

W = constantParam.mass*9.81;
rho = constantParam.rho;
R = constantParam.R_tip;
A_fus = constantParam.S_eq;

Vi_hover = sqrt(W/(2*rho*pi*(R^2)));

V = 1:0.5:100;

v_i = zeros(1, length(V));

for k = 1:length(V)
V_bar = V(k)/Vi_hover;
Drag_fus = 0.5 * rho * V(k)^2 * A_fus; %Uses equivalent flat plate area
alpha_d = asin(Drag_fus/W);





syms vi_D

eqn_vi_D = vi_D^4 +  (V_bar^2 + 2*V_bar*(sin(alpha_d)))*vi_D^2 - 1 == 0;
%eqn_vi_D = vi_D^4 +  (V_bar^2*vi_D^2) - 1 == 0;

vi_D = double(solve(eqn_vi_D,vi_D));

vi_D=vi_D(real(vi_D)>0&imag(vi_D)==0);

v_i(:, k) = vi_D*Vi_hover;

% syms vi_solve
% 
% eqn_vi_solve = vi_solve^4 + (vi_solve^2)*(2*V(k)*sin(alpha_d)) - 1 == 0;
% 
% solution = double(solve(eqn_vi_solve,vi_solve));
% 
% %solution = vpasolve(eqn_vi_solve, vi_solve);
% 
% vi_solution = max(solution(solution == real(solution)));
% 
% v_i(k) = vi_solution*Vi_hover;
end

figure
%plot(V, v_i, 'DisplayName', 'Induced velocity', LineWidth=1.0)
plot(V, v_i, LineWidth=1.0)
hold on
%legend('Interpreter','latex')
grid
title('Induced velocity in forward flight', 'Interpreter','latex')
xlabel('Forward velocity [m/s]')
xlim([0 100])
ylabel('Induced velocity [m/s]')

end