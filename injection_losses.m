function [deltaP_inj] = injection_losses(rho, m_dot, A_inj, C_d)
%INJECTION_LOSSES This function returns the pressure drop due to injection
%plates.

% PROTOTYPE:
%   [deltaP_inj] = injection_losses(rho, m_dot, A_inj, C_d)
%
% INPUT:
%  rho          Propellant density [kg/m^3]
%  m_dot        Propellant mass flow rate [kg/m^3]
%  A_inj        Area of one propellant injector [m^2]
%  Cd           Discharge coefficient [-]
%
% OUTPUT:
%  deltaP_inj  Pression drop due to injection plates [Pa]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-17:     First version

deltaP_inj = 8/rho*(m_dot/(A_inj*C_d))^2;

end

