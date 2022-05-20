function T = calculate_thrust(d_inj_ox, d_inj_f, A_th, Pc, Cd, rho_ox, rho_f, A_e)

Delta_P_inj = 0.2*Pc;

A_inj_ox = (pi*d_inj_ox^2)/4;
A_inj_f = (pi*d_inj_f^2)/4;
m_dot_ox = 2*A_inj_ox*Cd*sqrt(2*Delta_P_inj*rho_ox);
m_dot_f = A_inj_f*Cd*sqrt(2*Delta_P_inj*rho_f);

OF = m_dot_ox/m_dot_f;

p = CEA_int_C_F(A_e,A_th);

ct = Val_Int_C_F(p, Pc*10^(-5), OF);

T = A_th * Pc * ct;

end