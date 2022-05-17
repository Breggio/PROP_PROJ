clc
clear all 

rho_t = 1.1931; %[kg/m^3]
v_t = 1061.2; %[m/s]
D_t = 0.5692*1e-2 ; % [m]
r_t = D_t /2; %[m]
D_c  = 1.3905*1e-2; %[m]
r_c = D_c/2;  %[m]
mu_t = 0.000090209; % [Pa*s]

Re_t = (rho_t*v_t*D_t)/mu_t;

Re_mod = sqrt(r_t/r_c)*Re_t