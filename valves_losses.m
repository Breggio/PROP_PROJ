function [R_valves] = valves_losses(rho, A)
%VALVES_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to the valves (latch valve and check valve).
%
% PROTOTYPE:
%   [R_valves] = valves_losses(rho, A)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%  A            Cross-sectional wetted area [m^2]
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

% rho_H2O = 1000; % water density [kg/m^3]
% SG = rho/rho_H2O; % [-] specific gravity
% 
% %% Computation of C_v for each valve
% % Latch valve (open-close valve) p.9 on datasheet
% deltaP_latch = 19; % [psi] Pressure drop
% m_dot_latch = 0.4/0.14; % volumetric mass flow rate [gallons/min]
% % (lb/s -> gallons/min)
% 
% C_v_latch = m_dot_latch*sqrt(SG/deltaP_latch);
% 
% % Check valve
% deltaP_check = 3.0; % [psi] Pressure drop
% m_dot_check = 3.0/0.14; % [gallons/min] Volumetric mass flow rate
% 
% C_v_check = m_dot_check*sqrt(SG/deltaP_check);
% 
% %% Computation of K_v for each valve
% K_v_latch = 0.865*C_v_latch; % [m^3/h] Flow factor
% K_v_check = 0.865*C_v_check; % [m^3/h] Flow factor
% 
% % C_v_open_close = 13; % flow coefficient for open-close valve
% % C_v_check = 1.5; % flow coefficient for check valve
% 
% %% Computation of resistance coefficients
% R_latch = (3600^2/rho)*1/(rho_H2O*K_v_latch^2);
% R_check = (3600^2/rho)*1/(rho_H2O*K_v_check^2);
% 
% R_valves = R_latch + R_check;

%% Loss coefficients of valves
K_L_latch = 0.15;
K_L_check = 2;

%% Computation of resistance coefficient
R_latch = (0.5*K_L_latch)/(rho*A^2); % [Pa/(kg^2/s^2)]
R_check = (0.5*K_L_check)/(rho*A^2); % [Pa/(kg^2/s^2)]

R_valves = R_latch + R_check;

end

