function [R_cooling] = cooling_losses(f_D, d_h, L, rho, A)
%COOLING_LOSSES This function returns the coefficient that multiplies
% the square of the mass flow rate in the expression of the pressure drop
% due to regenerative cooling.
%
% PROTOTYPE:
%   [R_cooling] = cooling_losses(f_D, d_h, L, rho, A)
%
% INPUT:
%  f_D          Darcy friction factor [-]
%  d_h          Hydraulic diameter of the pipe [m]
%               (for a pipe of circular section, this equals D;
%               otherwise d_h = 4A/P for a pipe of cross-sectional area A
%               and perimeter P)          
%  L            Length of the pipe [m]
%  rho          Propellant density [kg/m^3]
%  A            Cross-sectional wetted area [m^2]
%
% OUTPUT:
%  R_cooling    Resistance coefficient for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_cooling = 0.5*f_D/d_h*L/(rho*A^2);

end

