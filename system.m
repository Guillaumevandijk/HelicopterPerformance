clear all 
close all


gamma = 9;
rho = 1.225;

Mass = 7000*9.81;
I_y = 100;
C_DtimesS = 1000;
Omega =  258    *2*pi/60;
R = 16.36/2;
sigma = 0.02;
Cl_alpha = 0.2;

%state variables:

u = 20;         %m/s
w = 3;          %m/s
q = 1;          %degrees/s
theta = 2;      %degrees

%control variables:
theta_0 = 4;    %degrees
theta_c = 3;    %degrees

%system
V = sqrt(u^2+w^2);
alpha_c = theta_c - atan(w/u);
lampda_c = V*sin(alpha_c)/(Omega*R);
mu = V*cos(alpha_c)/(Omega*R);


lampda_i = 0;   %choose
count = 0;
Flampda = 1; %dummy

while abs(Flampda)>0.001 && count<100

    lampda_i_2 = lampda_i+0.001;

    a_1 = (-16/gamma*(q*pi/180)/Omega  +8/3*mu*theta-2*mu*(lampda_i+lampda_c))/       (1-1/(2*mu^2));
    
    V_glau = V/(Omega*R)*cos(alpha_c-a_1);
    Vi_glau = V/(Omega*R)*sin(alpha_c-a_1)+lampda_i;
    
    CT_glau = 2*lampda_i*sqrt(V_glau^2+Vi_glau^2);
    
    CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i));
    
    Flampda = CT_elem - CT_glau;

    
    a_1 = (-16/gamma*(q*pi/180)/Omega  +8/3*mu*theta-2*mu*(lampda_i+lampda_c))/       (1-1/(2*mu^2));
    
    V_glau = V/(Omega*R)*cos(alpha_c-a_1);
    Vi_glau = V/(Omega*R)*sin(alpha_c-a_1)+lampda_i_2;
    
    CT_glau = 2*lampda_i*sqrt(V_glau^2+Vi_glau^2);
    
    CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i_2));
    
    Flampda_2 = CT_elem - CT_glau;

    Flampda_grad = (Flampda_2-Flampda)/(lampda_i_2-lampda_i);

   
    lampda_i = lampda_i - Flampda/Flampda_grad;


    count = count+1;
end








