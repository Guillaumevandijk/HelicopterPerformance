clear all 
close all

%Constants
%flight conditions
rho = 1.225;


%Helicopter parameters
mass =  8322.4;             %kg
chord = 0.5273;             %m
I_y = 87905.35;             %kgm^2  estimated
Omega =  257.83   *2*pi/60; %rad/s
R_tip = 8.1778;             %m
Cl_alpha = 0.1176   *180/pi;%1/rad
h = 2.4;                    %m, height from blades to cg
gamma = 9;                  %[-] still needs to be determined

%drag           still needs to be determined
Cd = 0;
S = 0;

%get solidity
A_blades = (R_tip)*chord*4;
A_disc = pi*R_tip^2;
sigma = A_blades/A_disc;


%% start parameters for simulation
%initial conditions
V_x0 =      46.3;  
theta_f0 =  -2   *pi/180;        %rad
q0 =        0   *pi/180;        %rad/s
u0 =        cos(theta_f0)*V_x0;   %m/s
w0 =        sin(theta_f0)*V_x0;   %m/s

%control variables:
theta_0 =   6   *pi/180;        %rad
theta_c =   6   *pi/180;        %rad

%tracking the heli start point
x = [0];
y = [0];

%timestep
dt = 0.01;

%% Create the system

%state variables:
u = u0;
w = w0;
theta_f = theta_f0;  
q = q0;


%starting simulation

sim_count = 0;
for i = 1:25

    
    
    %flight velocity
    V = sqrt(u^2+w^2);  
    %control plane angle
    alpha_c = theta_c - atan(w/u);
    alpha_c = theta_c + atan(w/u); %Lijkt logischer toch?
    
    %w fraction
    lampda_c = V*sin(alpha_c)/(Omega*R_tip);
    %u fraction
    mu = V*cos(alpha_c)/(Omega*R_tip);
    
    %% Find lampda_i, induced velocity fraction
    alpha_newton = 0.4; %newton raphson parameter
    
    %search start lampda
    lampda_i = 0;       %chosen
    count = 0;          %counting
    Flampda = 1;        %dummy to initiate while loop
    
    
    %trakcing the heli progression 
    V_x = cos(theta_f)*u+sin(theta_f)*w;
    V_y = cos(theta_f)*w-sin(theta_f)*u;
    x = [x x(end)+dt*V_x];
    y = [y y(end)+dt*V_y];

    
    while abs(Flampda)>0.0001 && count<1000
        lampda_i_2 = lampda_i+0.0001;
        
    
    
        %f(x)
        %obtain glauert thrust
        a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i+lampda_c))/       (1-1/2*mu^2);
        %absolute and induced velocity
        V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
        Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i;%glauert thrust
        CT_glau = 2*lampda_i*sqrt(V_glau^2+Vi_glau^2);%element thrust
        CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i));
        %match the two thrusts
        Flampda = CT_elem - CT_glau;
    
        %f(x+dx)
        %obtain glauert thrust
        a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i_2+lampda_c))/       (1-1/2*mu^2);
        %absolute and induced velocity
        V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
        Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i_2;
        CT_glau = 2*lampda_i_2*sqrt(V_glau^2+Vi_glau^2);%glauert thrust
        CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c - lampda_i_2));%element thrust
        %match the two thrusts
        Flampda_2 = CT_elem - CT_glau;
    
        %f'(x)
        Flampda_grad = (Flampda_2-Flampda)/(lampda_i_2-lampda_i);
    
        %xn+1 = x - alpha * f(x)/f'(x)
        lampda_i = lampda_i - alpha_newton* Flampda/Flampda_grad;
        count = count+1;
    end
    
    %get thrust and drag
    T = CT_glau*rho*(Omega*R_tip)^2*pi*R_tip^2;
    D = Cd*1/2*rho*V^2*S;
    
    %get derivatives
    u_dot = -9.81*sin(theta_f) - D*u/mass/V+T/mass*sin(theta_c-a_1) - q*w;
    w_dot = 9.81*cos(theta_f)- D*w/mass/V - T/mass*cos(theta_c-a_1) + q*u;
    q_dot = -T/I_y*h*sin(theta_c-a_1);
    theta_f_dot = q;
    
    %get state variables
    u = u+u_dot*dt;
    w = w+w_dot*dt;
    q = q+q_dot*dt;
    theta_f = theta_f+ theta_f_dot*dt;

    %count
    sim_count = sim_count+1;
    disp(sim_count);
end


%plot xy
plot(x,y);


% Try to figure out gamma
% Cl_alpha = 0.1176*180/pi;
% n = 200;
% dr  = 1/n; %Radial increment - adimensional
% r = (0:dr:1); %preallocating the r-range
% massMeter = 13.92; %Mass for meter kg/m
% massMeter = 243.94/R_tip;
% Ib = sum(massMeter*R_tip*((r*R_tip).^2)*dr); %Inertial Moment
% Ib = 243.94*(R_tip/2)^2;
% 
% % 
% gamma = rho*Cl_alpha*chord*R_tip^4/Ib; %doesnt hold up, why?
% 