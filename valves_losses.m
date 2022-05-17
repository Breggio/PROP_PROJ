function [R_valves] = valves_losses(rho)
%VALVES_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to the valves (latch valve and check valve).
%
% PROTOTYPE:
%   [R_valves] = valves_losses(rho)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%
% OUTPUT:
%  R_valves     Resistance coefficient for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

rho_H2O = 1000; % water density [kg/m^3]

C_v_open_close = 13; % flow coefficient for open-close valve
C_v_check = 1.5; % flow coefficient for check valve

K_v_open_close = 0.865*C_v_open_close; % flow factor for open-close valve [m^3/h]
K_v_check = 0.865*C_v_check; % flow factor for check valve [m^3/h]

R_open_close = (3600^2/rho)*1/(rho_H2O*K_v_open_close^2);
R_check = (3600^2/rho)*1/(rho_H2O*K_v_check^2);

R_valves = R_open_close + R_check;

end

