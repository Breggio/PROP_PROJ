function [c_star] = Val_Int_C_star(p, P_c, OF)
%VAL_INT_C_STAR This function returns the value of the characteristic
% velocity evaluating the interpolated polynomial for given combustion
% chamber pressure and OF ratio.
%
% PROTOTYPE:
%   [c_star] = Val_Int_C_star(p, P_c, OF)
%
% INPUT:
%  p            Coefficients of the polynomial obtained by interpolation
%  P_c          Combustion chamber pressure [bar]
%  OF           Oxidizer to fuel mass ratio [-]
%
% OUTPUT:
%  c_star       Characteristic velocity [m/s]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-20: First version

c_star = polyVal2D(p,P_c,OF);

end

