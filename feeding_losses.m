function [deltaP_feed] = feeding_losses(f, rho, L, D, m_dot)
%FEEDING_LOSSES This function returns the pressure drop due to feeding
%line.
%
% PROTOTYPE:
%   [deltaP_feed] = feeding_losses(f, rho, L, D, m_dot)
%
% INPUT:
%  f            Darcy friction factor [-]
%  rho          Propellant density [kg/m^3]
%  L            Length of the feeding line [m]
%  D            Diameter of the pipe in the feeding line [m]
%  m_dot        Propellant mass flow rate [kg/m^3]
%
% OUTPUT:
%  deltaP_feed  Pression drop due to feeding line [Pa]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

deltaP_feed = 8/pi^2*f*m_dot^2/rho*L/D^5; % [Pa] Pressure drop

end

