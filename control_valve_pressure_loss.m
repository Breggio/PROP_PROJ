function pressure_loss = control_valve_pressure_loss(volumic_mass)
%Returns the pressure loss of a control valve at our conditions, given the
%loss given in the datasheet of the valve.
%   CAREFUL: GIVE NUMBERS IN THE DEMANDED UNITS (else results WILL be wrong)
%   INPUTS : 
%   volumic_mass : the volumic mass of the fluid used in the real
%   system (kg/m³)
%   OUTPUT : Rlatch : factor multiplying mdot² to give DeltaP
%   (bar.kg^-2.s^2)
.
%   VERSIONS : 2022/05/19, Alexandre BOULEUX
%
Rlatch = 1.31 * (1/(0.0001814))^2 * (1/volumic_mass*1000);
Rcheck = 0.21 * (1/(0.0012))^2 * (1/volumic_mass*1141);
pressure_loss = Rlatch+Rcheck
end