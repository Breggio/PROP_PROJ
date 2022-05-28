close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Create log(Rb) vs log(P) plot

% Load pressure data and calculate Rb values

load tracesbar1.mat;

plow_1 = pbar2438(:,1);
pmid_1 = pbar2438(1:4201,2);
phigh_1 = pbar2438(1:3601,3);

Rb = zeros(9,3);
Peff = zeros(9,3);

[Rb(1,1), Peff(1,1), tburn1(1,1), tburn2(1,1)] = calulate_Rb(plow_1);
[Rb(1,2), Peff(1,2), tburn1(1,1), tburn2(1,1)] = calulate_Rb(pmid_1);
[Rb(1,3), Peff(1,3), tburn1(1,1), tburn2(1,1)] = calulate_Rb(phigh_1);

plow_2 = pbar2439(:,1);
pmid_2 = pbar2439(1:4201,2);
phigh_2 = pbar2439(1:3601,3);

[Rb(2,1), Peff(2,1), tburn1(2,1), tburn2(2,1)] = calulate_Rb(plow_2);
[Rb(2,2), Peff(2,2), tburn1(2,2), tburn2(2,2)] = calulate_Rb(pmid_2);
[Rb(2,3), Peff(2,3), tburn1(2,3), tburn2(2,3)] = calulate_Rb(phigh_2);

plow_3 = pbar2440(:,1);
pmid_3 = pbar2440(1:4201,2);
phigh_3 = pbar2440(1:3601,3);

[Rb(3,1), Peff(3,1), tburn1(3,1), tburn2(3,1)] = calulate_Rb(plow_3);
[Rb(3,2), Peff(3,2), tburn1(3,2), tburn2(3,2)] = calulate_Rb(pmid_3);
[Rb(3,3), Peff(3,3), tburn1(3,3), tburn2(3,3)] = calulate_Rb(phigh_3);

plow_4 = pbar2441(:,1);
pmid_4 = pbar2441(1:4301,2);
phigh_4 = pbar2441(1:3601,3);

[Rb(4,1), Peff(4,1), tburn1(4,1), tburn2(4,1)] = calulate_Rb(plow_4);
[Rb(4,2), Peff(4,2), tburn1(4,2), tburn2(4,2)] = calulate_Rb(pmid_4);
[Rb(4,3), Peff(4,3), tburn1(4,3), tburn2(4,3)] = calulate_Rb(phigh_4);

plow_5 = pbar2442(:,1);
pmid_5 = pbar2442(1:4201,2);
phigh_5 = pbar2442(1:3601,3);

[Rb(5,1), Peff(5,1), tburn1(5,1), tburn2(5,1)] = calulate_Rb(plow_5);
[Rb(5,2), Peff(5,2), tburn1(5,2), tburn2(5,2)]= calulate_Rb(pmid_5);
[Rb(5,3), Peff(5,3), tburn1(5,3), tburn2(5,3)] = calulate_Rb(phigh_5);

plow_6 = pbar2443(:,1);
pmid_6 = pbar2443(1:4301,2);
phigh_6 = pbar2443(1:3601,3);

[Rb(6,1), Peff(6,1), tburn1(6,1), tburn2(6,1)] = calulate_Rb(plow_6);
[Rb(6,2), Peff(6,2), tburn1(6,2), tburn2(6,2)] = calulate_Rb(pmid_6);
[Rb(6,3), Peff(6,3), tburn1(6,3), tburn2(6,3)] = calulate_Rb(phigh_6);

plow_7 = pbar2444(:,1);
phigh_7 = pbar2444(1:4601,2); % QUESTI DUE SONO INVERTITIIIIIIIIIIIIIIIIIIII
pmid_7 = pbar2444(1:4201,3);

pmid_7(10) = 5;

figure()
plot([1:length(pmid_7)], pmid_7);

[Rb(7,1), Peff(7,1), tburn1(7,1), tburn2(7,1)] = calulate_Rb(plow_7);
[Rb(7,2), Peff(7,2), tburn1(7,2), tburn2(7,2)] = calulate_Rb(pmid_7);
[Rb(7,3), Peff(7,3), tburn1(7,3), tburn2(7,3)] = calulate_Rb(phigh_7);

plow_8 = pbar2445(:,1);
pmid_8 = pbar2445(1:4201,2);
phigh_8 = pbar2445(1:3601,3);

[Rb(8,1), Peff(8,1), tburn1(8,1), tburn2(8,1)] = calulate_Rb(plow_8);
[Rb(8,2), Peff(8,2), tburn1(8,2), tburn2(8,2)] = calulate_Rb(pmid_8);
[Rb(8,3), Peff(8,3), tburn1(8,3), tburn2(8,3)] = calulate_Rb(phigh_8);

plow_9 = pbar2446(:,1);
pmid_9 = pbar2446(1:4301,2);
phigh_9 = pbar2446(1:3601,3);

[Rb(9,1), Peff(9,1), tburn1(9,1), tburn2(9,1)] = calulate_Rb(plow_9);
[Rb(9,2), Peff(9,2), tburn1(9,2), tburn2(9,2)] = calulate_Rb(pmid_9);
[Rb(9,3), Peff(9,3), tburn1(9,3), tburn2(9,3)] = calulate_Rb(phigh_9);

Peff_vec = reshape(Peff, 27, 1);
Rb_vec = reshape(Rb, 27, 1);

c = polyfit(log(Peff_vec), log(Rb_vec), 1);
Rb_line = polyval(c, log(Peff_vec));

figure()
for k=1:9
    for j=1:3
        plot(log(Peff(k,j)), log(Rb(k,j)), 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
        hold on
    end
end
plot(log(Peff_vec), Rb_line, 'k', 'LineWidth', 2.5);
grid on
grid minor
xlabel('log(P)')
ylabel('log($r_{b}$)')

n = c(1); % [-]
a = exp(c(2)); % [mm/(s bar^n)]
a_Pa = a/(1e5)^n; % [m/(s Pa^n)]

[a_test, Inc_a, n_test, Inc_n, R2] = Uncertainty(Peff_vec, Rb_vec);
% Inc_a [mm/(s bar^n)]

% Throat area of the engines
At_low = 28.8*1e-3; % [m]
At_mid = 25.25*1e-3; % [m]
At_high = 21.81*1e-3; % [m]

% Propellant density
% ERANO SBAGLIATE LE UNITà DI MISURA. M_ox/f/b erano [kg] non [g]

M_ox = 68; % [kg]
rho_ox = 1950; % [g/cm^3]
M_f = 18; % [kg]
rho_f = 2700; % [g/cm^3]
M_b = 14; % [kg]
rho_b = 920; % [g/cm^3]
rho_p = (M_ox + M_f + M_b)/(M_ox/rho_ox + M_f/rho_f + M_b/rho_b); % [kg/m^3]

% Characteristic velocity (amount of mass exiting from a combustion chamber)

% Ideal
MM = 25.574; % [g/mol]
Cp = 2.0456*MM; % [KJ/(kg*K)]
R = 8.314472; % [KJ/(mol*K)]
gamma = Cp/(Cp-R);
m_mol = MM*10^(-3);
R_specific = R/m_mol;
T = 3333.74; % [K]
c_star = sqrt(gamma*R_specific*T)*1/(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1))));

% % Real
% V_p = pi*((0.08)^2-(0.05)^2)*0.29; % [m^3]
% Mtot = rho_p*V_p; % [kg]
% At_matrix = zeros(9,3);
% At_matrix(1,1) = At_low;
% At_matrix(:,2) = At_mid;
% At_matrix(:,3) = At_high;
% clear i
% clear j
% for i=1:9
%     for j=1:3
%         c_star_matrix(i,j) = Peff(i,j)*(tburn2(i,j) - tburn1(i,j))*1e5*At_matrix(i,j)/Mtot;
%     end
% end
% c_star_low = mean(c_star_matrix(:,1));
% c_star_mid = mean(c_star_matrix(:,2));
% c_star_high = mean(c_star_matrix(:,3));

% Burning rate
Rb_low = mean(Rb(:,1));
Rb_mid = mean(Rb(:,2));
Rb_high = mean(Rb(:,3));

% Set the time step
delta_t = 0.001; % [sec]

% Set the time vector
Time = 80*60*1e3; % [sec]

c_star_low = 1610.4;
c_star_mid = 1616.6;
c_star_high = 1622.8;

[tb_low, P_low] = BARIA(a_test, n, delta_t, At_low, c_star_low, rho_p, Time);
[tb_mid, P_mid] = BARIA(a_test, n, delta_t, At_mid, c_star_mid, rho_p, Time);
[tb_high, P_high] = BARIA(a_test, n, delta_t, At_high, c_star_high, rho_p, Time);

%% Monte Carlo

nb_samples = 50;
nb_iterations = 10;

for it = 1:nb_iterations

% 1) Define the population
a_vec = random('norm',a,Inc_a,nb_samples,1);
n_vec = random('norm',n,Inc_n,nb_samples,1);

% 2) Forming the couples
Couples = zeros(1, nb_samples*nb_samples);

z = 1;
for i = 1:nb_samples
    for j = 1:nb_samples
            Couples(z) = a_vec(i);
            Couples(z+1) = n_vec(j);
            z = z+2;
    end
end

% 4) Shuffling couples
nb_couples = length(Couples)/2;

indexes_couples = (1:nb_couples); % creating a list with indexes of couples
indexes_couples = indexes_couples(randperm(nb_couples));
% randomly shuffle the indexes of the couples

Couples_shuffled = zeros(1, length(Couples));

for i = 1:nb_couples
    couples_index = indexes_couples(i);
    Couples_shuffled(2*i-1) = Couples(2*couples_index-1);
    Couples_shuffled(2*i) = Couples(2*couples_index);
end

% 5) Run the simulation for each couple

tb_low_vec = zeros(1, nb_samples*nb_samples);
tb_mid_vec = zeros(1, nb_samples*nb_samples);
tb_high_vec = zeros(1, nb_samples*nb_samples);

k = 1;
for i = 1:length(Couples_shuffled)
    if mod(i,2) == 1 % selecting only odd indexes QUESTO NON SO CHE CAZZO FA
        couple = Couples_shuffled(i:i+1);
        a_mc = couple(1);
        n_mc = couple(2);
        [tb_low, P_low] = BARIA(a_mc, n_mc, delta_t, At_low, c_star_low, rho_p, Time);
        [tb_mid, P_mid] = BARIA(a_mc, n_mc, delta_t, At_mid, c_star_mid, rho_p, Time);
        [tb_high, P_high] = BARIA(a_mc, n_mc, delta_t, At_high, c_star_high, rho_p, Time);
        tb_low_vec(k) = tb_low;
        tb_mid_vec(k) = tb_mid;
        tb_high_vec(k) = tb_high;
        k = k+1;
    end
end

mean_tb_low = mean(tb_low_vec);
mean_tb_mid = mean(tb_mid_vec);
mean_tb_high = mean(tb_high_vec);

standard_deviation_tb_low = std(tb_low_vec);
standard_deviation_tb_mid = std(tb_mid_vec);
standard_deviation_tb_high = std(tb_high_vec);

Mean_tb_low(it) = mean_tb_low;
Mean_tb_mid(it) = mean_tb_mid;
Mean_tb_high(it) = mean_tb_high;

Standard_deviation_tb_low(it) = standard_deviation_tb_low;
Standard_deviation_tb_mid(it) = standard_deviation_tb_mid;
Standard_deviation_tb_high(it) = standard_deviation_tb_high;

Cumulative_mean_tb_low(it) = mean(Mean_tb_low(1:it));
Cumulative_mean_tb_mid(it) = mean(Mean_tb_mid(1:it));
Cumulative_mean_tb_high(it) = mean(Mean_tb_high(1:it));

Cumulative_standard_deviation_tb_low(it) = mean(Standard_deviation_tb_low(1:it));
Cumulative_standard_deviation_tb_mid(it) = mean(Standard_deviation_tb_mid(1:it));
Cumulative_standard_deviation_tb_high(it) = mean(Standard_deviation_tb_high(1:it));

end

% 6) Map the code convergence
figure();

% Plotting cumulative mean
subplot(1,2,1);
plot(Cumulative_mean_tb_low, 'LineWidth', 1.5);
hold on
plot(Cumulative_mean_tb_mid, 'LineWidth', 1.5);
plot(Cumulative_mean_tb_high, 'LineWidth', 1.5);
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative mean [s]', 'Interpreter', 'latex');
title('\textbf{Cumulative mean of burning time}', 'Interpreter', 'latex');
legend('Low P', 'Mid P', 'High P')

% Plotting cumulative standard deviation
subplot(1,2,2);
plot(Cumulative_standard_deviation_tb_low, 'LineWidth', 1.5);
hold on
plot(Cumulative_standard_deviation_tb_mid, 'LineWidth', 1.5);
plot(Cumulative_standard_deviation_tb_high, 'LineWidth', 1.5);
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative standard deviation [s]', 'Interpreter', 'latex');
title('\textbf{Cumulative standard deviation of burning time}', 'Interpreter', 'latex');
legend('Low P', 'Mid P', 'High P')

sgtitle('\textbf{Monte Carlo analysis}', 'Interpreter', 'latex');

%% 7) Determine the criterion for convergence of cumulative mean and std deviation (95% coverage interval)

% Sort the Monte Carlo results from lowest value to the highest value
Cumulative_mean_tb_low_sorted = sort(Cumulative_mean_tb_low);
Cumulative_mean_tb_mid_sorted = sort(Cumulative_mean_tb_mid);
Cumulative_mean_tb_high_sorted = sort(Cumulative_mean_tb_high);

Cumulative_standard_deviation_tb_low_sorted = sort(Cumulative_standard_deviation_tb_low);
Cumulative_standard_deviation_tb_mid_sorted = sort(Cumulative_standard_deviation_tb_mid);
Cumulative_standard_deviation_tb_high_sorted = sort(Cumulative_standard_deviation_tb_high);

% For a 95% coverage interval
if 0.025*nb_iterations == fix(0.025*nb_iterations) % if 0.025*nb_iterations is an integer
    cm_low_low = Cumulative_mean_tb_low_sorted(0.025*nb_iterations);
    cm_low_mid = Cumulative_mean_tb_mid_sorted(0.025*nb_iterations);
    cm_low_high = Cumulative_mean_tb_high_sorted(0.025*nb_iterations);

    sd_low_low = Cumulative_standard_deviation_tb_low_sorted(0.025*nb_iterations);
    sd_low_mid = Cumulative_standard_deviation_tb_mid_sorted(0.025*nb_iterations);
    sd_low_high = Cumulative_standard_deviation_tb_high_sorted(0.025*nb_iterations);
else
    index = fix(0.025*nb_iterations + 0.5);
    cm_low_low = Cumulative_mean_tb_low_sorted(index);
    cm_low_mid = Cumulative_mean_tb_mid_sorted(index);
    cm_low_high = Cumulative_mean_tb_high_sorted(index);

    sd_low_low = Cumulative_standard_deviation_tb_low_sorted(index);
    sd_low_mid = Cumulative_standard_deviation_tb_mid_sorted(index);
    sd_low_high = Cumulative_standard_deviation_tb_high_sorted(index);
end

if 0.975*nb_iterations == fix(0.975*nb_iterations)  % if 0.975*nb_iterations is an integer
    cm_high_low = Cumulative_mean_tb_low_sorted(0.975*nb_iterations);
    cm_high_mid = Cumulative_mean_tb_mid_sorted(0.975*nb_iterations);
    cm_high_high = Cumulative_mean_tb_high_sorted(0.975*nb_iterations);

    sd_high_low = Cumulative_standard_deviation_tb_low_sorted(0.975*nb_iterations);
    sd_high_mid = Cumulative_standard_deviation_tb_mid_sorted(0.975*nb_iterations);
    sd_high_high = Cumulative_standard_deviation_tb_high_sorted(0.975*nb_iterations);
else
    index = fix(0.975*nb_iterations + 0.5);
    cm_high_low = Cumulative_mean_tb_low_sorted(index);
    cm_high_mid = Cumulative_mean_tb_mid_sorted(index);
    cm_high_high = Cumulative_mean_tb_high_sorted(index);

    sd_high_low = Cumulative_standard_deviation_tb_low_sorted(index);
    sd_high_mid = Cumulative_standard_deviation_tb_mid_sorted(index);
    sd_high_high = Cumulative_standard_deviation_tb_high_sorted(index);
end

% For 95% expanded limit uncertainty
mean_cm_low = mean(Cumulative_mean_tb_low_sorted);
mean_cm_mid = mean(Cumulative_mean_tb_mid_sorted);
mean_cm_high = mean(Cumulative_mean_tb_high_sorted);

mean_sd_low = mean(Cumulative_standard_deviation_tb_low_sorted);
mean_sd_mid = mean(Cumulative_standard_deviation_tb_mid_sorted);
mean_sd_high = mean(Cumulative_standard_deviation_tb_high_sorted);

U_cm_minus_low = mean_cm_low - cm_low_low;
U_cm_minus_mid = mean_cm_mid - cm_low_mid;
U_cm_minus_high = mean_cm_high - cm_low_high;

U_cm_plus_low = cm_high_low - mean_cm_low;
U_cm_plus_mid = cm_high_mid - mean_cm_mid;
U_cm_plus_high = cm_high_high - mean_cm_high;

U_sd_minus_low = mean_sd_low - sd_low_low;
U_sd_minus_mid = mean_sd_mid - sd_low_mid;
U_sd_minus_high = mean_sd_high - sd_low_high;

U_sd_plus_low = sd_high_low - mean_sd_low;
U_sd_plus_mid = sd_high_mid - mean_sd_mid;
U_sd_plus_high = sd_high_high - mean_sd_high;

% The interval that contains r_true at a 95% level of confidence is then
Int_cm_inf_low = mean_cm_low - U_cm_minus_low;
Int_cm_inf_mid = mean_cm_mid - U_cm_minus_mid;
Int_cm_inf_high = mean_cm_high - U_cm_minus_high;

Int_cm_sup_low = mean_cm_low + U_cm_plus_low;
Int_cm_sup_mid = mean_cm_mid + U_cm_plus_mid;
Int_cm_sup_high = mean_cm_high + U_cm_plus_high;

Int_sd_inf_low = mean_sd_low - U_sd_minus_low;
Int_sd_inf_mid = mean_sd_mid - U_sd_minus_mid;
Int_sd_inf_high = mean_sd_high - U_sd_minus_high;

Int_sd_sup_low = mean_sd_low + U_sd_plus_low;
Int_sd_sup_mid = mean_sd_mid + U_sd_plus_mid;
Int_sd_sup_high = mean_sd_high + U_sd_plus_high;

% Adding the convergence intervals to the plots
subplot(1,2,1);
yline([Int_cm_inf_low Int_cm_sup_low], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');
hold on
yline([Int_cm_inf_mid Int_cm_sup_mid], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');
yline([Int_cm_inf_high Int_cm_sup_high], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');

subplot(1,2,2);
yline([Int_sd_inf_low Int_sd_sup_low], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');
hold on
yline([Int_sd_inf_mid Int_sd_sup_mid], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');
yline([Int_sd_inf_high Int_sd_sup_high], '--', {'Lower bound','Upper bound'}, 'LineWidth', 1, 'color', '#D95319');
