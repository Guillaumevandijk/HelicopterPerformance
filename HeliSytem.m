clear all 
close all




[constantParam] = getConstants();
%Trimmed_condition(constantParam);

%Constants
%flight conditions
rho = constantParam.rho; %kg/m^2
%Helicopter parameters
mass = constantParam.mass;%kg
chord= constantParam.chord;            %m
I_y = constantParam.I_y;             %kgm^2  estimated
Omega = constantParam.Omega; %rad/s
R_tip = constantParam.R_tip;             %m
Cl_alpha = constantParam.Cl_alpha;%1/rad
h= constantParam.h;                    %m, height from blades to cg
gamma = constantParam.gamma;                  %[-] still needs to be determined
%drag           still needs to be determined
Cd = constantParam.Cd;
S = constantParam.S;
%get solidity
A_blades = constantParam.A_blades;
A_disc = constantParam.A_disc;
sigma = constantParam.sigma;
A_fus =constantParam.S_eq;


%% start parameters for simulation
%initial conditions
V_x0 =      0.01; 
V_y0 =      0;

%get theta_f0 with drag
D = 0.5 * rho * V_x0^2 * A_fus;
theta_f0 = atan(-D/(mass*9.81))   *pi/180;        %rad
q0 =        0   *pi/180;        %rad/s
u0 =        cos(theta_f0)*V_x0;   %m/s
w0 =        sin(theta_f0)*V_x0;   %m/s

%control variables:
theta_0 =   8.94   *pi/180;        %rad
theta_c =  1   *pi/180;        %rad

%tracking the heli start point
x_track = 0;
y_track = 100;
V_x_track = V_x0;
V_y_track = V_y0;
u_track = u0;
w_track = w0;
theta_f_track = theta_f0;
time_track = 0;

%timestep
dt = 0.001;

%% Create the system

%state variables:
u = u0;
w = w0;
theta_f = theta_f0;  
q = q0;


%starting simulation

time = 0;
disp(" alright");
check = true;
sum = 0;
while time <3 
    

    %flight velocity
    V = sqrt(u^2+w^2);  
    
    %control plane angle
    if u >= 0
        alpha_c = theta_c - theta_f;
    else
        alpha_c = theta_c + theta_f;
    end

    %w fraction
    lampda_c = V*sin(alpha_c)/(Omega*R_tip);
    %u fraction
    mu = V*cos(alpha_c)/(Omega*R_tip);
    
    %% Find lampda_i, induced velocity fraction
    alpha_newton = 0.05; %newton raphson parameter
    
    %search start lampda
    lampda_i = 0;       %chosen
    count = 0;          %counting
    Flampda = 1;        %dummy to initiate while loop
   
    

    while abs(Flampda)>0.0001 && count<2000
        lampda_i_2 = lampda_i+0.0000001;
        
        if count>1999
        disp(time);
        
        end
    
        %f(x)
        %obtain glauert thrust
        a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i+lampda_c))/       (1-1/2*mu^2);
        %absolute and induced velocity
        V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
        Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i;%glauert thrust
        CT_glau = 2*lampda_i*sqrt(V_glau^2+Vi_glau^2);%element thrust
        CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c + lampda_i));
        %match the two thrusts
        Flampda = CT_elem - CT_glau;
    
        %f(x+dx)
        %obtain glauert thrust
        a_1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i_2+lampda_c))/       (1-1/2*mu^2);
        %absolute and induced velocity
        V_glau = V/(Omega*R_tip)*cos(alpha_c-a_1);
        Vi_glau = V/(Omega*R_tip)*sin(alpha_c-a_1)+lampda_i_2;
        CT_glau = 2*lampda_i_2*sqrt(V_glau^2+Vi_glau^2);%glauert thrust
        CT_elem = 1/4*Cl_alpha*sigma*(2/3*theta_0*(1+3/2*mu^2)-(lampda_c + lampda_i_2));%element thrust
        %match the two thrusts
        Flampda_2 = CT_elem - CT_glau;
    
        %f'(x)
        Flampda_grad = (Flampda_2-Flampda)/(lampda_i_2-lampda_i);
    
        %xn+1 = x - alpha * f(x)/f'(x)
        lampda_i = lampda_i - alpha_newton* Flampda/2/Flampda_grad;
        count = count +1;
    end
    
    %get thrust and drag
    T = CT_glau*rho*(Omega*R_tip)^2*pi*R_tip^2;
    D = 0.5 * rho * V^2 * A_fus; %Uses equivalent flat plate area
    
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
    time = time+dt;
    %disp(time);

    %tracking the heli progression 
    V_x = cos(-theta_f)*u-sin(-theta_f)*w;
    V_y = -cos(-theta_f)*w-sin(-theta_f)*u;

    x_track = [x_track x_track(end)+dt*V_x];
    y_track = [y_track y_track(end)+dt*V_y];
    V_x_track = [V_x_track V_x];
    V_y_track = [V_y_track V_y];
    time_track = [time_track time];
    u_track = [u_track u];
    w_track = [w_track w];
    theta_f_track = [theta_f_track theta_f];

    %theta_c = theta_c -  0.000002*abs((V_x-30));

end

%plot xy

% Create a figure and subplot layout
figure(1);
subplot(4,1,1); 
plot(-x_track, y_track);
title('coordinates');
xlabel('-X');
ylabel('Height');
axis([-200,200,0,400]);


subplot(4,1,2); 
plot(time_track,V_x_track );
title('Vx over time');
xlabel('t');
ylabel('Vx');

subplot(4,1,3); 
plot(time_track,V_y_track);
title('Vy over time');
xlabel('t');
ylabel('Vy');

subplot(4,1,4); 
plot(time_track,180/pi*theta_f_track);
title('theta_f over time');
xlabel('t');
ylabel('theta_f');




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