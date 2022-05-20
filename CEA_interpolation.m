function [outputArg1,outputArg2] = CEA_interpolation(b, a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A_e = 21.4179e-3; % [m^2]
A_t = 2.6772e-5; % [m^2] 

load('CEAmatrix.mat');
load('Pe.mat');

%% Replacing 4th column of v with Pe values
v(:,4) = Pe';

%% Constants
R = 8.314; % molar gas constant [kJ/(K.kmol)]
nb_inputs = 5;

OF_vect = linspace(6.8, 8, nb_inputs);
Pc_vect = linspace(20, 10, nb_inputs);

P_c = [Pc_vect; Pc_vect; Pc_vect; Pc_vect; Pc_vect]; % c.c. pressure [bar]
T = zeros(nb_inputs, nb_inputs); % matrix of temperatures [K]
M_mol = zeros(nb_inputs, nb_inputs); % matrix of molar masses [kg/kmol]
C_p = zeros(nb_inputs, nb_inputs); % matrix of specific heat capacities [kJ/(kg.K)]
C_star = zeros(nb_inputs, nb_inputs); % c* [m/s]
C_F = zeros(nb_inputs, nb_inputs); % thrust coefficient [-]
P_e = zeros(nb_inputs, nb_inputs); % exit pressure [bar]

% T: 2nd column, M_mol: 3rd, P_e: 4th, C_p: 5th
for i = 1:size(v, 2)
    for j = 1:size(v, 2)
        index = 5*(i-1)+j;
        T(i,j) = v(index,2);
        M_mol(i,j) = v(index,3);
        P_e(i,j) = v(index, 4);
        C_p(i,j) = v(index,5)*M_mol(i,j);
    end
end

Gamma = C_p./(C_p-R); % [-]

%% C_star and C_F computation
eps = A_e/A_t;

for i = 1:size(C_star, 1)
    for j = 1:size(C_star, 2)
        %C_star
        t = T(i,j);
        m_mol = M_mol(i,j)*10^(-3); % molar mass [kg/mol]
        R_specific = R/m_mol;  % specific gas constant [J/(kg.K)]
        gamma = Gamma(i,j);
        c_star = sqrt(gamma*R_specific*t)*1/(gamma*sqrt((2/(gamma+1))^((gamma+1)/(gamma-1))));
        C_star(i,j) = c_star;

        %C_F
        k = gamma;
        Pe2Pc = P_e(i,j)/P_c(i,j);
        c_F = sqrt(2*(k^2/(k-1))*(2/(k+1))^((k+1)/(k-1)))*sqrt(1-(Pe2Pc)^((k-1)/k))+eps*Pe2Pc;
        C_F(i,j) = c_F;
    end
end



end