function pressure_loss = control_valve_pressure_loss(given_pressure_loss, given_flow_rate, give_volumic_mass, mdot, real_volumic_mass)
%Returns the pressure loss of a control valve at our conditions, given the
%loss given in the datasheet of the valve.
%   CAREFUL: GIVE NUMBERS IN THE DEMANDED UNITS (else results WILL be wrong)
%   INPUTS : 
%   given_pressure_loss : the pressure drop given in the datasheet
%   (bar, or else problems)
%   given_flow_rate : the volumetric flow rate for the given pressure drop
%   in the datasheet (m³/s)
%   given_volumic_mass : the volumic mass of the fluid used in the
%   datasheet (usually water, but check) (kg/m³)
%   mdot : real mass flow rate (kg/s)
%   real_volumic_mass : the volumic mass of the fluid used in the real
%   system (kg/m³)
%   OUTPUT : pressure_loss : the pressure drop induced by the valve (bar)
%   VERSIONS : 2022/05/19, Alexandre BOULEUX
%
real_flow_rate = mdot/real_volumic_mass
pressure_loss = given_pressure_loss * (real_flow_rate/given_flow_rate)^2 * (real_volumic_mass/give_volumic_mass);

end