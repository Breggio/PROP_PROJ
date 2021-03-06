close all
clear
clc

%% Constants
P_c = 20E5; % nominal combustion chamber pressure [Pa]
C_d = 0.7; % nominal discharge coefficient [-]
%OF = 7.2; % nominal oxidizer to fuel ratio [-]
rho_f = 810; % fuel density [kg/m^3]
rho_ox = 1373; % oxidizer density [kg/m^3]
A_e = 0.002071; % nozzle exit area [m^2]

%m_dot = 0.0327; % propellant mass flow rate [kg/s]
%m_dot_ox = (OF/(1+OF))*m_dot; % nominal oxidizer mass flow rate [kg/s]
%m_dot_f = m_dot - m_dot_ox; % nominal fuel mass flow rate [kg/s]

%% Preparing vectors for the results
nb_iterations = 10; % number of iterations

Mean_thrust = zeros(1, nb_iterations);
Cumulative_mean_thrust = zeros(1, nb_iterations);
Standard_deviation_thrust = zeros(1, nb_iterations);
Cumulative_standard_deviation_thrust = zeros(1, nb_iterations);

%% 1) Parameters
nb_samples = 5; % number of samples

nom_d_inj_ox = 8.8721e-04; % nominal value of oxidizer injector diameter [m]
nom_d_inj_f = 5.3355e-04; % nominal value of fuel injector diameter [m]
tol_d_inj = 7.5e-05; % tolerance on injector diameter [m]
nom_A_th = 2.5884e-05; % nominal value of throat area [m^2]
tol_A_th = (pi*(9E-5)^2)/4; % tolerance on throat area [m^2]

%% STARTING LOOP OF ITERATIONS
for it = 1:nb_iterations

%% 2) Define the population
% Injection hole diameter oxidizer
a_d_inj_ox = nom_d_inj_ox - tol_d_inj;
b_d_inj_ox = nom_d_inj_ox + tol_d_inj;
D_ox = (b_d_inj_ox-a_d_inj_ox).*rand(1,nb_samples) + a_d_inj_ox;

% Injection hole diameter fuel
a_d_inj_f = nom_d_inj_f - tol_d_inj;
b_d_inj_f = nom_d_inj_f + tol_d_inj;
D_f = (b_d_inj_f-a_d_inj_f).*rand(1,nb_samples) + a_d_inj_f;

% Throat diameter
a_A_th = nom_A_th - tol_A_th;
b_A_th = nom_A_th + tol_A_th;
A_th = a_A_th + (b_A_th-a_A_th).*rand(1,nb_samples);

%% 3) Forming the couples
Triplets = zeros(1, nb_samples*nb_samples*nb_samples);

z = 1;
for i = 1:nb_samples
    for j = 1:nb_samples
        for k = 1:nb_samples
            Triplets(z) = D_ox(i);
            Triplets(z+1) = D_f(j);
            Triplets(z+2) = A_th(k);
            z = z+3;
        end
    end
end

%% 4) Shuffling couples
nb_triples = length(Triplets)/3;

indexes_triplets = (1:nb_triples); % creating a list with indexes of triplets
indexes_triplets = indexes_triplets(randperm(nb_triples));
% randomly shuffle the indexes of the triplets

Triplets_shuffled = zeros(1, length(Triplets));

for i = 1:nb_triples
    triplet_index = indexes_triplets(i);
    Triplets_shuffled(3*i-2) = Triplets(3*triplet_index-2);
    Triplets_shuffled(3*i-1) = Triplets(3*triplet_index-1);
    Triplets_shuffled(3*i) = Triplets(3*triplet_index);
end

%% 5) Run the simulation for each triplet
Thrust = zeros(1, nb_samples*nb_samples*nb_samples);

k = 1;
for i = 1:length(Triplets_shuffled)
    if mod(i,3) == 1 % selecting only odd indexes
        triplet = Triplets_shuffled(i:i+2);
        d_ox = triplet(1);
        d_f = triplet(2);
        a_th = triplet(3);
        thrust = calculate_thrust(d_ox, d_f, a_th, P_c, C_d, rho_ox, rho_f, A_e);
        Thrust(k) = thrust;
        k = k+1;
    end
end

mean_thrust = mean(Thrust);
standard_deviation_thrust = std(Thrust);

Mean_thrust(it) = mean_thrust;
Standard_deviation_thrust(it) = standard_deviation_thrust;

Cumulative_mean_thrust(it) = mean(Mean_thrust(1:it));
Cumulative_standard_deviation_thrust(it) = mean(Standard_deviation_thrust(1:it));

end
%% CLOSING LOOP OF ITERATIONS

%% 6) Map the code convergence
figure(1);

% Plotting cumulative mean
subplot(1,2,1);
plot(Cumulative_mean_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative mean [N]', 'Interpreter', 'latex');
title('\textbf{Cumulative mean of thurst}', 'Interpreter', 'latex');

% Plotting cumulative standard deviation
subplot(1,2,2);
plot(Cumulative_standard_deviation_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative standard deviation [N]', 'Interpreter', 'latex');
title('\textbf{Cumulative standard deviation of thrust}', 'Interpreter', 'latex');

sgtitle('\textbf{Monte Carlo analysis}', 'Interpreter', 'latex');

%% 7) Determine the criterion for convergence of cumulative mean and std deviation (95% coverage interval)

% Sort the Monte Carlo results from lowest value to the highest value
Cumulative_mean_thrust_sorted = sort(Cumulative_mean_thrust);
Cumulative_standard_deviation_thrust_sorted = sort(Cumulative_standard_deviation_thrust);

% For a 95% coverage interval
if 0.025*nb_iterations == fix(0.025*nb_iterations) % if 0.025*nb_iterations is an integer
    cm_low = Cumulative_mean_thrust_sorted(0.025*nb_iterations);
    sd_low = Cumulative_standard_deviation_thrust_sorted(0.025*nb_iterations);
else
    index = fix(0.025*nb_iterations + 0.5);
    cm_low = Cumulative_mean_thrust_sorted(index);
    sd_low = Cumulative_standard_deviation_thrust_sorted(index);
end

if 0.975*nb_iterations == fix(0.975*nb_iterations)  % if 0.975*nb_iterations is an integer
    cm_high = Cumulative_mean_thrust_sorted(0.975*nb_iterations);
    sd_high = Cumulative_standard_deviation_thrust_sorted(0.975*nb_iterations);
else
    index = fix(0.975*nb_iterations + 0.5);
    cm_high = Cumulative_mean_thrust_sorted(index);
    sd_high = Cumulative_standard_deviation_thrust_sorted(index);
end

% For 95% expanded limit uncertainty
mean_cm = mean(Cumulative_mean_thrust_sorted);
mean_sd = mean(Cumulative_standard_deviation_thrust_sorted);

U_cm_minus = mean_cm - cm_low;
U_cm_plus = cm_high - mean_cm;

U_sd_minus = mean_sd - sd_low;
U_sd_plus = sd_high - mean_sd;

% The interval that contains r_true at a 95% level of confidence is then
Int_cm_inf = mean_cm - U_cm_minus;
Int_cm_sup = mean_cm + U_cm_plus;

Int_sd_inf = mean_sd - U_sd_minus;
Int_sd_sup = mean_sd + U_sd_plus;

% Adding the convergence intervals to the plots
subplot(1,2,1);
yline([Int_cm_inf Int_cm_sup], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');

subplot(1,2,2);
yline([Int_sd_inf Int_sd_sup], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');

