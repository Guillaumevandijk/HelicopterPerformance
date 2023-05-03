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

%drag           still needs to be determined
constantParam.Cd = 0;
constantParam.S = 0;

%get solidity
constantParam.A_blades = (constantParam.R_tip)*constantParam.chord*4;
constantParam.A_disc = pi*constantParam.R_tip^2;
constantParam.sigma = constantParam.A_blades/constantParam.A_disc;
