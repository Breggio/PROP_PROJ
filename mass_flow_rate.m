function [m_dot] = mass_flow_rate(P_gt, P_c, R_tot)
%MASS_FLOW_RATE This function returns the mass flow rate taking into
%account the losses through the feeding system, the pressure of the gas
%inside the tank and the combustion chamber pressure.
%
% PROTOTYPE:
%   [m_dot] = mass_flow_rate(P_gt, P_c, R_tot)
%
% INPUT:
%  P_gt         Gas pressure in tank [Pa]
%  P_c          Combustion chamber pressure [Pa]
%  R_tot        Sum of all resistance coefficients for propellant flow path
%               [Pa/(kg^2/s^2)]
%
% OUTPUT:
%  m_dot        Mass flow rate [kg/s]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

m_dot = sqrt( (P_gt-P_c)/R_tot ); % [kg/s]

end

