function [R_feed] = feeding_losses(f, rho, L, D)
%FEEDING_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to feeding line.
%
% PROTOTYPE:
%   [R_feed] = feeding_losses(f, rho, L, D)
%
% INPUT:
%  f            Darcy friction factor [-]
%  rho          Propellant density [kg/m^3]
%  L            Length of the feeding line [m]
%  D            Diameter of the pipe in the feeding line [m]
%
% OUTPUT:
%  R_feed       Resistance coefficient for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_feed = 8/pi^2*f/rho*L/D^5; % [Pa/(kg^2/s^2)]

end

