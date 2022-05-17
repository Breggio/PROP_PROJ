function [m_dot] = mass_flow_rate(P_gt, P_c, R_feed, R_dyn, R_inj, R_valve)
%MASS_FLOW_RATE This function returns the mass flow rate taking into
%account the losses through the feeding system, the pressure of the tank
%and the combustion chanmber pressure.
%
% PROTOTYPE:
%   [m_dot] = mass_flow_rate(P_gt, P_c, R_feed, R_dyn, R_inj, R_valve)
%
% INPUT:
%  P_gt         Pressure
%  P_c          Combustion chamber pressure [Pa]
%  R_feed       Pressure drop coefficient due to feeding line [Pa/(kg^2/s^2)]
%  R_dyn        Pressure drop coefficient due to dynamic losses [Pa/(kg^2/s^2)]
%  R_inj        Pressure drop coefficient due to injection plate [Pa/(kg^2/s^2)]
%  R_valve      Pressure drop coefficient due to valve [Pa/(kg^2/s^2)]
%
% OUTPUT:
%  m_dot        Mass flow rate [kg/s]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

R_sum = R_feed + R_dyn + R_inj + R_valve; % [Pa/(kg^2/s^2)]
m_dot = sqrt( (P_gt-P_c)/ R_sum); % [kg/s]

end

