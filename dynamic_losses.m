function [deltaP_dyn] = dynamic_losses(rho, m_dot, A)
%DYNAMIC_LOSSES This function returns the pressure drop due at the exit of
%the tank.
%
% PROTOTYPE:
%   [deltaP_dyn] = dynamic_losses(rho, m_dot, A)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%  m_dot        Propellant mass flow rate [kg/m^3]
%  A            Area of the pipe at the exit of the tank [m^2]
%
% OUTPUT:
%  deltaP_dyn   Pression drop at the exit of the tank [Pa]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

deltaP_dyn = 1/2*m_dot^2/(rho*A^2);

end

