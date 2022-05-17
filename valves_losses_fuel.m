function [R_valves] = valves_losses_fuel(rho)
%VALVES_LOSSES Summary of this function goes here
%   Detailed explanation goes here
rho_H2O = 1000; % [kg/m^3]

C_v_open = 13; % flow coefficient for open-close valve
C_v_check = 1.5; % flow coefficient for check valve

K_v_open = 0.865*C_v_open; % flow factor for open-close valve [m^3/h]
K_v_check = 0.865*C_v_check; % flow factor for check valve [m^3/h]

R_open = (3600^2/rho)/(rho_H2O*K_v_open^2);
R_check = (3600^2/rho)/(rho_H2O*K_v_check^2);

R_valves = R_open + R_check;

end