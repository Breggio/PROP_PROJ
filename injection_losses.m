function [R_inj] = injection_losses(rho, A_inj, C_d)
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
%
% OUTPUT:
%  R_inj        Pression drop coefficient [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_inj = (8/rho)/(A_inj*C_d)^2;

end

