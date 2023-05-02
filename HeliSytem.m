clear all 
close all



rho = 1.225;

Mass =  8322.4*9.81;
chord = 0.5273;
I_y = 13693;
C_DtimesS = 1000;
Omega =  257.83   *2*pi/60;
sigma = 0.0681;
Cl_alpha = 0.1176*180/pi;


R_tip = 8.1778;
A_blades = (R_tip)*chord*4;
A_disc = pi*R_tip^2;
sigma = A_blades/A_disc;

n = 200;
dr  = 1/n; %Radial increment - adimensional
r = (0:dr:1); %preallocating the r-range
massMeter = 13.92; %Mass for meter kg/m
massMeter = 243.94/R_tip;
Ib = sum(massMeter*R_tip*((r*R_tip).^2)*dr); %Inertial Moment
Ib = 243.94*(R_tip/2)^2;

% 
gamma = rho*Cl_alpha*chord*R_tip^4/Ib; %doesnt hold up, why?

gamma = 9;




theta_f = 2*pi/180;      %rad
V_x = 46.3;

%state variables:
u = cos(theta_f)*V_x; %m/s
w = sin(theta_f)*V_x; %m/s
q = 0   *pi/180;          %rad/s


%control variables:
theta_0 = 6*pi/180;    %rad
theta_c = 6*pi/180;    %rad

%system
V = sqrt(u^2+w^2);
alpha_c = theta_c - atan(w/u);
lampda_c = V*sin(alpha_c)/(Omega*R_tip)
mu = V*cos(alpha_c)/(Omega*R_tip);

lampda_i = 0;   %choose
count = 0;
Flampda = 1; %dummy
alpha_newton = 0.2; %newton raphson parameter
while abs(Flampda)>0.001 && count<1000

    lampda_i_2 = lampda_i+0.001;

    a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i+lampda_c))/       (1-1/2*mu^2);
    V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
    Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i;
    CT_glau = 2*lampda_i*sqrt(V_glau^2+Vi_glau^2);
    CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i));
    Flampda = CT_elem - CT_glau;

    
    a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i_2+lampda_c))/       (1-1/2*mu^2);
    V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
    Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i_2;
    CT_glau = 2*lampda_i_2*sqrt(V_glau^2+Vi_glau^2);
    CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i_2));
    Flampda_2 = CT_elem - CT_glau;



    Flampda_grad = (Flampda_2-Flampda)/(lampda_i_2-lampda_i);

   
    lampda_i = lampda_i - alpha_newton* Flampda/Flampda_grad;


    count = count+1;
end








