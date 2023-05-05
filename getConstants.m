function [constantParam] = getConstants()
%Constants
%flight conditions
constantParam.rho = 1.225;


%Helicopter parameters
constantParam.mass =  8322.4;             %kg
constantParam.chord = 0.5273;             %m
constantParam.I_y = 87905.35;             %kgm^2  estimated
constantParam.Omega =  257.83   *2*pi/60; %rad/s
constantParam.R_tip = 8.1778;             %m
constantParam.Cl_alpha = 0.1176   *180/pi;%1/rad
constantParam.h = 2.4;                    %m, height from blades to cg
constantParam.gamma = 9;                  %[-] still needs to be determined
constantParam.FM = 0.8;


%drag           still needs to be determined
constantParam.Cd = 0;
constantParam.S = 0;
constantParam.S_eq = 3.1773; %[m^2]
constantParam.C_dp_main = 0.0057;

%get solidity
constantParam.A_blades = (constantParam.R_tip)*constantParam.chord*4;
constantParam.A_disc = pi*constantParam.R_tip^2;
constantParam.sigma = constantParam.A_blades/constantParam.A_disc;

constantParam.k_main = 1.1;
constantParam.k_tail = 1.3;
constantParam.Tailarm = 9; %[m] schatting
constantParam.R_tail = 1.68; %[m]


%NOG INVULLEN
constantParam.C_dp_tail = 0.0057;
constantParam.chord_tail = 0.12;             %m
constantParam.Omega_tail = 257.83   *2*pi/60; %rad/s 
constantParam.A_blades_tail = (constantParam.R_tail)*constantParam.chord_tail*4;
constantParam.A_disc_tail = pi*constantParam.R_tail^2;
constantParam.sigma_tail = constantParam.A_blades_tail/constantParam.A_disc_tail;


