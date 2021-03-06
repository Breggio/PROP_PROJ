clc
clear all

rho_t = 1.1931; %[kg/m^3]
v_t = 1061.2; %[m/s]
D_t = 0.5692*1e-2 ; % [m]
r_t = D_t /2; %[m]
D_c  = 1.3905*1e-2; %[m]
r_c = D_c/2;  %[m]
mu_t = 0.000090209; % [Pa*s]
mu_c = 0.000096284 % [Pa*s]
Re_t = (rho_t*v_t*D_t)/mu_t;

Re_mod = sqrt(r_t/r_c)*Re_t;

%% feeding line 

D_pipe_ox = 10e-3; %m
D_pipe_fu = 5e-3; %m
rho_fu = 810; %kg/m3
mu_fu = 1.40e-3; %pa*s


rho_ox = 1373; %kg/m3 
mu_ox = 1.26e-3; %pa*s 
v_fu = 0.2507;  %m/s
v_ox = 0.2663;  %m/s

Re_ox = (rho_ox * D_pipe_ox * v_ox)/mu_ox;
Re_fu = (rho_fu * D_pipe_fu * v_fu)/mu_fu; 