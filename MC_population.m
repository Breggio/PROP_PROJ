close all
clear
clc

%% 1) Parameters
N_samples = 2; % number of samples

nom_d_inj_ox = 9.2114e-04; % nominal value of oxidizer injector diameter [m]
nom_d_inj_f = 5.5395e-04; % nominal value of fuel injector diameter [m]
tol_d_inj = 0.08e-03; % tolerance on injector diameter [m]
nom_A_th = 2.5450e-05; % nominal value of throat area [m^2]
tol_A_th = (10e-06)^2; % tolerance on throat area [m^2]

%% 2) Define the population
% Injection hole diameter oxidizer
a_d_inj_ox = nom_d_inj_ox - tol_d_inj;
b_d_inj_ox = nom_d_inj_ox + tol_d_inj;
D_ox = (b_d_inj_ox-a_d_inj_ox).*rand(1,N_samples) + a_d_inj_ox;

% Injection hole diameter fuel
a_d_inj_f = nom_d_inj_f - tol_d_inj;
b_d_inj_f = nom_d_inj_f + tol_d_inj;
D_f = (b_d_inj_f-a_d_inj_f).*rand(1,N_samples) + a_d_inj_f;

% Throat diameter
a_A_th = nom_A_th - tol_A_th;
b_A_th = nom_A_th + tol_A_th;
A_th = a_A_th + (b_A_th-a_A_th).*rand(1,N_samples);

%% 3) Forming the couples
Triplets = zeros(1, N_samples*N_samples*N_samples);

z = 1;
for i = 1:N_samples
    for j = 1:N_samples
        for k = 1:N_samples
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
Thrust = zeros(1, N_samples*N_samples*N_samples);

k = 1;
for i = 1:length(Triplets_shuffled)
    if mod(i,3) == 1 % selecting only odd indexes
        triplet = Triplets_shuffled(i:i+1);

        Thrust(k) = thrust;
        k = k+1;
    end
end