function T = calculate_thrust(d_inj_ox, d_inj_f, A_th, Pc, Cd, rho_ox, rho_f)

Delta_P_inj= 0.25*Pc;

A_inj_ox = (pi*d_inj_ox^2)/4;
A_inj_f = (pi*d_inj_f^2)/4;
m_dot_ox = 2*A_inj_ox*Cd*sqrt(2*Delta_P_inj*rho_ox);
m_dot_f = A_inj_f*Cd*sqrt(2*Delta_P_inj*rho_f);

OF = m_dot_ox/m_dot_f;

[outputs] = CEA('problem','rocket','frozen','o/f',OF,'case','CEAM-rocket1',...
    'p,Pa',Pc,'supsonic(ae/at)',80,'reactants','fuel','RP-1(L)','C',1,...
    'H',1.95000,'wt%',100,'t(k)',298.0,'oxid','H2O2(L)','wt%',87.5,...
    't(k)',350,'oxid','H2O(L)','wt%',12.5,'t(k)',350,...
    'output','thermochemical','end');

ct = outputs.output.froz.cf_vac(3);

T = A_th * Pc * ct;

end