close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Thrust Monte Carlo Simulation

%% Constants

%% Variable ranges

% % Roughness combustion chamber
% N_f_cc = 10;
% f_cc_vec = a + (b-a).*rand(N_f_cc,1);
% QUESTO NEL CASO CAPIAMO COME CALCOLARE LA PERDITA DI PRESSIONE DOVUTA
% ALLA ROUGHNESS SUPERFICIALE IN CAMERA


% Injection holes diameter
N_d_inj = 10;
d_inj_vec = 0.5E-3  (b-a).*rand(N_inj_vec,1);
% Throat diameter
N_d_th = 10;
d_th_vec = a + (b-a).*rand(N_th_vec,1);
% Combustion chamber diameter
N_d_cc = 10;
d_cc_vec = a + (b-a).*rand(N_d_cc,1);


