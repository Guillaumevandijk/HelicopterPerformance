clear all;
close all;

R = 16.36/2;
W = 7000*9.81;
alpha_c = 2     *pi/180;
Omega =  258    *2*pi/60;
gamma = 9;
theta = 2       *pi/180;
rho = 1.225;

q = 20          *pi/180;
V = 300;

T = W/cos(alpha_c);

vi_h = sqrt(T/(2*rho*R^2));
V_D = V/vi_h;



syms vi_D

eqn_vi_D = vi_D^4 +V_D^2*vi_D^2-1 == 0;

vi_D = double(solve(eqn_vi_D,vi_D));

 vi_D=vi_D(real(vi_D)>0&imag(vi_D)==0);

v_i = vi_D*vi_h;



lampda_i = v_i/(Omega*R);
lampda_c = V*sin(alpha_c)/(Omega*R);
mu = V*cos(alpha_c)/(Omega*R);


a0 = (gamma/8)*     (theta* (1+mu^2-  4/3*(lampda_i+lampda_c) )   );
a1 = (-16/gamma*(q)/Omega  +8/3*mu*theta-2*mu*(lampda_i+lampda_c))/       (1-1/(2*mu^2));
b1 = -q/Omega+4/3*mu*a0/    (1+1/(2*mu^2));


%% plot blades

circP = 30;
radP = 40;

psi = [0:pi/circP:2*pi]';
beta = a0 - a1*cos(psi)-b1*sin(psi);
r = 0:R/radP:R;
z = tan(beta)*r;
x = sin(psi)*r;
y = -cos(psi)*r;



xlabel("X");
ylabel("Y");
zlabel("Z");
%surf(x,y,z)


%% 


beta_dot = Omega*(a1*sin(psi) - b1*cos(psi));
alpha = theta - v_i./r/Omega- beta_dot/Omega + q/Omega*cos(psi);


y = sin(psi)*r  %x, y switched for view
x = -cos(psi)*r;

xlabel("Xdsafadsf");
ylabel("Y");
zlabel("Z");
contour(x,y,alpha,23);


% view(0,90)
% %make hub plane
% r = 0:R/radP/2:R/2;
% beta = alpha_c*cos(psi)-b1*sin(psi);
% z = tan(beta)*r;
% x = sin(psi)*r;
% y = -cos(psi)*r;

