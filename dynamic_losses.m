function [R_dyn] = dynamic_losses(rho, A)
%DYNAMIC_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to dynamic losses (at the exit of the tank).
%
% PROTOTYPE:
%   [R_dyn] = dynamic_losses(rho, A)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%  A            Area of the pipe at the exit of the tank [m^2]
%
% OUTPUT:
%  R_dyn        Resistance coefficient for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_dyn = 0.5/(rho*A^2);

end

