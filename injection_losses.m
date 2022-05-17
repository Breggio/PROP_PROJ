function [R_inj] = injection_losses(rho, A_inj, C_d, N_inj)
%INJECTION_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to injection plate.

% PROTOTYPE:
%   [R_inj] = injection_losses(rho, A_inj, C_d)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%  A_inj        Area of one propellant injector [m^2]
%  Cd           Discharge coefficient [-]
%  N_inj        Number of injectors
%
% OUTPUT:
%  R_inj        Resistance coefficient for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_inj = N_inj*(8/rho)/(A_inj*C_d)^2;

% if R_inj > 0.3*P_c
%     disp("Attention: the injection losses are too high!")
% else
%     disp("Injection losses are within 5-30% of P_c! BRAVO!!")
% end

end

