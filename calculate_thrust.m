function T = calculate_thrust(d_inj_ox, d_inj_f, A_th, Pc_old, OF, Cd, rho_ox, rho_f, m_dot_ox, m_dot_f)

Delta_P_inj_old = 0.2*Pc_old;

A_inj_ox = (pi*d_inj_ox^2)/4;
Delta_P_inj_ox_new = (1/(2*rho_ox))*(m_dot_ox/(A_inj_ox*Cd));

A_inj_f = (pi*d_inj_f^2)/4;
Delta_P_inj_f_new = (1/(2*rho_f))*(m_dot_f/(A_inj_f*Cd));

Pc_new = Pc_old - (Delta_P_inj_f_new - Delta_P_inj_old) - (Delta_P_inj_ox_new - Delta_P_inj_old);

[outputs] = CEA('problem','rocket','frozen','o/f',OF,'case','CEAM-rocket1',...
    'p,Pa',Pc_new,'supsonic(ae/at)',80,'reactants','fuel','RP-1(L)','C',1,...
    'H',1.95000,'wt%',100,'t(k)',298.0,'oxid','H2O2(L)','wt%',87.5,...
    't(k)',350,'oxid','H2O(L)','wt%',12.5,'t(k)',350,...
    'output','thermochemical','end');

ct = outputs.output.froz.cf_vac(3);

T = A_th * Pc_new * ct;

end