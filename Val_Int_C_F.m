function [c_F] = Val_Int_C_F(p, P_c, OF)
%VAL_INT_C_STAR This function returns the value of the thrust coefficient
% evaluating the interpolated polynomial for given combustion chamber
% pressure and OF ratio.
%
% PROTOTYPE:
%   [c_F] = Val_Int_C_F(p, P_c, OF)
%
% INPUT:
%  p            Coefficients of the polynomial obtained by interpolation
%  P_c          Combustion chamber pressure [bar]
%  OF           Oxidizer to fuel mass ratio [-]
%
% OUTPUT:
%  c_F          Thrust coefficient [-]
%
% CONTRIBUTORS:
%   Léonie DEU
%
% VERSIONS
%   2022-05-20: First version

c_F = polyVal2D(p,P_c,OF,2,2);

end