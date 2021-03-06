close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Constants
P_c = 20E5; % nominal combustion chamber pressure [Pa]
C_d = 0.65; % nominal discharge coefficient [-]
OF = 7.2; % nominal oxidizer to fuel ratio [-]
rho_f = 810; % fuel density [kg/m^3]
rho_ox = 1373; % oxidizer density [kg/m^3]

m_dot = 0.0327; % propellant mass flow rate [kg/s]
m_dot_ox = (OF/(1 + OF))*m_dot; % [kg/s] Oxidizer mass flow rate
m_dot_f = m_dot - m_dot_ox; % [kg/s] Fuel mass flow rate

%% 1) Parameters
nb_samples = 10; % number of samples

nom_d_inj_ox = 9.2114e-04; % nominal value of oxidizer injector diameter [m]
nom_d_inj_f = 5.5395e-04; % nominal value of fuel injector diameter [m]
tol_d_inj = 0.08e-03; % tolerance on injector diameter [m]
nom_A_th = 2.5450e-05; % nominal value of throat area [m^2]
tol_A_th = (10e-06)^2; % tolerance on throat area [m^2]

%% 2) Define the population
% Injection hole diameter oxidizer
a_d_inj_ox = nom_d_inj_ox - tol_d_inj;
b_d_inj_ox = nom_d_inj_ox + tol_d_inj;
D_ox_vec = (b_d_inj_ox-a_d_inj_ox).*rand(nb_samples,1) + a_d_inj_ox;

% Injection hole diameter fuel
a_d_inj_f = nom_d_inj_f - tol_d_inj;
b_d_inj_f = nom_d_inj_f + tol_d_inj;
D_f_vec = (b_d_inj_f-a_d_inj_f).*rand(nb_samples,1) + a_d_inj_f;

% Throat diameter
a_A_th = nom_A_th - tol_A_th;
b_A_th = nom_A_th + tol_A_th;
A_th_vec = a_A_th + (b_A_th-a_A_th).*rand(nb_samples,1);

%% 3) INitialize the thrust matrix
[T_array] = deal(zeros(nb_samples, nb_samples, nb_samples));

%% 4) MC run

% if no visibility problems and Parallel Computing Toolbox present,
% outermost for can be replaced by "parfor" to speed up
for i = 1 : nb_samples
    for j = 1 : nb_samples
        for k = 1 : nb_samples
            d_inj_ox = D_ox_vec(i);
            d_inj_f = D_f_vec(j);
            A_th = A_th_vec(k);
            T = calculate_thrust(d_inj_ox, d_inj_f, A_th, P_c, C_d, rho_ox, rho_f);
            T_array(i,j,k) = T;
        end
    end
end

% Save the resulting matrix
save('MC_results.mat')

% Plot the 3D matrix

% d_inj_ox vector for the plot
d_inj_ox_plot = [D_ox_vec(1).*ones(100,1); D_ox_vec(2).*ones(100,1); D_ox_vec(3).*ones(100,1); D_ox_vec(4).*ones(100,1); D_ox_vec(5).*ones(100,1); D_ox_vec(6).*ones(100,1); D_ox_vec(7).*ones(100,1); D_ox_vec(8).*ones(100,1); D_ox_vec(9).*ones(100,1); D_ox_vec(10).*ones(100,1)];

% d_inj_f vector for the plot
d_inj_f_plot = repmat([D_f_vec(1).*ones(10,1); D_f_vec(2).*ones(10,1); D_f_vec(3).*ones(10,1); D_f_vec(4).*ones(10,1); D_f_vec(5).*ones(10,1); D_f_vec(6).*ones(10,1); D_f_vec(7).*ones(10,1); D_f_vec(8).*ones(10,1); D_f_vec(9).*ones(10,1); D_f_vec(10).*ones(10,1)],10,1);

% A_th vector for the plot
A_th_plot = repmat(A_th_vec,100,1);

% thrust vector for the plot
T_vec = reshape(T_array,nb_samples*nb_samples*nb_samples,1);

%%

% Plot
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

s = scatter3(d_inj_ox_plot, d_inj_f_plot, A_th_plot,100,T_vec,'filled');
a = colorbar
set(a, 'TickLabelInterpreter', 'latex');
xlabel('Ox injector diameter [m]')
ylabel('Fuel injector diameter [m]')
zlabel('Throat area [$m^2$]')
colorTitleHandle = get(a,'Title');
titleString = 'Thrust [N]';
set(colorTitleHandle ,'String',titleString,'Interpreter','latex');
colormap turbo

%%

% Maximum value
thrust_max = max(max(max(T_array)))

% Minimum value 
thrust_min = min(min(min(T_array)))

% Mean
thrust_mean = mean(T_array,'all')

% Standard deviation
thrust_std = std(T_array,0,'all')



