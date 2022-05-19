function pressure_loss = control_valve_pressure_loss(mdot, volumic_mass)
%Returns the pressure loss of a control valve at our conditions, given the
%loss given in the datasheet of the valve.
%   CAREFUL: GIVE NUMBERS IN THE DEMANDED UNITS (else results WILL be wrong)
%   INPUTS : 
%   mdot : real mass flow rate (kg/s)
%   real_volumic_mass : the volumic mass of the fluid used in the real
%   system (kg/m³)
%   OUTPUT : pressure_loss : factor multiplying mdot² to give DeltaP
%   (bar.kg^-2.s^2)
.
%   VERSIONS : 2022/05/19, Alexandre BOULEUX
%
pressure_loss = 1.31 * (1/(0.0001814*real_volumic_mass))^2 * (volumic_mass/1000);

end