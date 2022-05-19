close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Thrust Monte Carlo Simulation

%% Constants

%% Variable ranges

<<<<<<< Updated upstream
% % Roughness combustion chamber
% N_f_cc = 10;
% f_cc_vec = a + (b-a).*rand(N_f_cc,1);
% QUESTO NEL CASO CAPIAMO COME CALCOLARE LA PERDITA DI PRESSIONE DOVUTA
% ALLA ROUGHNESS SUPERFICIALE IN CAMERA


% Injection holes diameter
N_d_inj = 10;
d_inj_vec = 0.5E-3  (b-a).*rand(N_inj_vec,1);
=======
% Injection hole diameter oxidizer
N_d_inj_ox = 100;
a_d_inj_ox = (9.2114e-04)-(0.08e-03);
b_d_inj_ox = (9.2114e-04)+(0.08e-03);
d_inj_ox_vec = a_d_inj_ox + (b_d_inj_ox-a_d_inj_ox).*rand(N_inj_ox_vec,1);
% Injection hole diameter fuel
N_d_inj_f = 100;
a_d_inj_f = (5.5395e-04)-(0.08e-03);
b_d_inj_f = (5.5395e-04)+(0.08e-03);
d_inj_f_vec = a_d_inj_f + (b_d_inj_f-a_d_inj_f).*rand(N_inj_f_vec,1);
>>>>>>> Stashed changes
% Throat diameter
N_A_th = 100;
a_A_th = (2.5450e-05)-(10e-06);
b_A_th = (2.5450e-05)+(10e-06);
A_th_vec = a_A_th + (b_A_th-a_A_th).*rand(N_A_th,1);



