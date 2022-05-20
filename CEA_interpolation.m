function [outputArg1,outputArg2] = CEA_interpolation(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load('CEAmatrix.mat');

%% Constants
R = 8.314; % molar gas constant [J/(K.mol)]
nb_inputs = 5;

OF_vect = linspace(6.8, 8, nb_inputs);
Pc_vect = linspace(20, 10, nb_inputs);

T = zeros(nb_inputs, nb_inputs); % matrix of temperatures
M_mol = zeros(nb_inputs, nb_inputs); % matrix of molar masses
C_p = zeros(nb_inputs, nb_inputs); % matrix of specific heat capacities

% T: 2nd column, M_mol: 3rd, C_p: 5th
for i = 1:size(v, 1)
    for j = 1:size(v, 2)
        T(i,j) = v(j,2);
        M_mol(i,j) = v(j,3);
        C_p(i,j) = v(j,5);
    end
end




end