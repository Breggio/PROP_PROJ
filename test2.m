clear
clc

A_e = 2.14179e-3; % [m^2]
A_t = 2.6772e-5; % [m^2] 

load('CEAmatrix.mat');
load('Pe.mat');

%% Replacing 4th column of v with Pe values
v(:,4) = Pe';

%% Constants
R = 8.314; % molar gas constant [kJ/(K.kmol)]
nb_inputs = 7;
OF_low = 7;
OF_high = 8.3;
Pc_low = 20;
Pc_high = 5;

OF_vect = linspace(OF_low, OF_high, nb_inputs);
Pc_vect = linspace(Pc_low, Pc_high, nb_inputs);

%% Creating matrices to store the results
P_c = repmat(Pc_vect, nb_inputs); % c.c. pressure [bar]
T = zeros(nb_inputs, nb_inputs); % matrix of temperatures [K]
M_mol = zeros(nb_inputs, nb_inputs); % matrix of molar masses [kg/kmol]
C_p = zeros(nb_inputs, nb_inputs); % matrix of specific heat capacities [kJ/(kg.K)]
C_star = zeros(nb_inputs, nb_inputs); % c* [m/s]
C_F = zeros(nb_inputs, nb_inputs); % thrust coefficient [-]
P_e = zeros(nb_inputs, nb_inputs); % exit pressure [bar]

%% Filling matrices
% T: 2nd column, M_mol: 3rd, P_e: 4th, C_p: 5th
for i = 1:nb_inputs
    for j = 1:nb_inputs
        index = nb_inputs*(i-1)+j;
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

%% Interpolation
[X,Y] = meshgrid(Pc_vect, OF_vect);

p = polyFit2D(C_F,X,Y,2,2);

c_F = polyVal2D(p,20,7.8,2,2);
