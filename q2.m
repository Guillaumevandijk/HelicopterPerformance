clear all;
close all;

R = 
Omega = 
v_i = 
R = 
alpha_c = 
Omega =   *2*pi/60;
gamma = 
theta = 2;


q = 20

lampda_i = v_i/(Omega*R);
lampda_c = V*sind(alpha_c)/(Omega*R);
mu = V*cosd(alpha_c)/(Omega*R);


a0 = (gamma/8)*     (theta* (1+mu^2-  4/3*(lampda_i+lampda_c) )   );
a1 = (-16/gamma*(q*pi/180)/Omega  +8/3*mu*theta-2*mu*(lampda_i+lampda_c))/       (1-1/(2*mu^2));
b1 = -q/Omega+4/3*mu*a0/    (1+1/(2*mu^2));


beta = a0 - a1*cosd(psi)-b1*sind(psi);
