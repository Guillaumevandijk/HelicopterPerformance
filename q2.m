clear all;
close all;

%Obtain constants
constantParam = getConstants();
R_tip = constantParam.R_tip;
mass = constantParam.mass;
Omega =  constantParam.Omega;
gamma = constantParam.gamma;
rho = constantParam.rho;

%Variables
q = 20          *pi/180;
p = 10          *pi/180;
V = 20;
theta_0 = 6       *pi/180;
alpha_c = 2     *pi/180;
alpha_cl= 1     *pi/180;


%Hoover induced velocity
T = mass*9.81/cos(alpha_c);
vi_h = sqrt(T/(2*rho*R_tip^2));
V_D = V/vi_h;

%find V_i
syms vi_D
eqn_vi_D = vi_D^4 +V_D^2*vi_D^2-1 == 0;
vi_Dall = double(solve(eqn_vi_D,vi_D));
vi_D=vi_Dall(real(vi_Dall)>0&imag(vi_Dall)==0);
v_i = vi_D*vi_h;

%get velocity fractions
lampda_i = v_i/(Omega*R_tip);
lampda_c = V*sin(alpha_c)/(Omega*R_tip);
mu = V*cos(alpha_c)/(Omega*R_tip);

%find disc angles
a0 = (gamma/8)*     (theta_0* (1+mu^2)-  4/3*(lampda_i+lampda_c)   );
a0_s = (gamma/8)*     (theta_0-  4/3*lampda_i);

a1 = (-16/gamma*q/Omega  +8/3*mu*theta_0-2*mu*(lampda_i+lampda_c))/  (1-0.5*mu^2);
a1_s = -16/gamma*q/Omega;


b1 = (-q/Omega+4/3*mu*a0)/ (1+0.5*mu^2);
b1_s = -q/Omega;


%% plot blades

circP = 17;
radP = 15;


psi = linspace(0,2*pi,circP)';
beta = a0 - a1*cos(psi)-b1*sin(psi);
r = linspace(1.5,R_tip,radP);
z = tan(beta)*r;
x = sin(psi)*r;
y = -cos(psi)*r;

%surf(x,y,z)


%plot the AoA plot

beta_dot = Omega*(a1*sin(psi) - b1*cos(psi));
teller = V*sin(alpha_c)   +v_i   +  beta_dot.*r   -   q*r.*cos(psi)   +   V*cos(alpha_c).*cos(psi).*sin(beta);
noemer = Omega*r+V*cos(alpha_c)*sin(psi);

alpha = theta_0 - teller./noemer;


figure();

x = sin(psi)*r;  %x, y switched for view
y = -cos(psi)*r;

xlabel("Xdsafadsf");
ylabel("Y");
zlabel("Z");
axis([-9,9,-9,9,-inf,inf])

contour(x,y,alpha*180/pi,-10:1:4,'ShowText','on');
figure()
surf(x,y,alpha*180/pi);
view(0,45);























% view(0,90)
% %make hub plane
% r = 0:R_tip/radP/2:R_tip/2;
% beta = alpha_c*cos(psi)-b1*sin(psi);
% z = tan(beta)*r;
% x = sin(psi)*r;
% y = -cos(psi)*r;
% 
% xlabel("X");
% ylabel("Y");
% zlabel("Z");
% 
% %surf(x,y,z);



%% 

