


d_inj_ox_low = a_d_inj_ox ;
d_inj_f_low = a_d_inj_f;

Delta_P_inj_old = 0.2*Pc;

A_inj_ox = (pi*d_inj_ox^2)/4;
Delta_P_inj_ox_new = (1/(2*rho_ox))*(m_dot_ox/(A_inj_ox*Cd));

A_inj_f = (pi*d_inj_f^2)/4;
Delta_P_inj_f_new = (1/(2*rho_f))*(m_dot_f/(A_inj_f*Cd));

Pc_new_low = Pc_old - (Delta_P_inj_f_new - Delta_P_inj_old) - (Delta_P_inj_f_new - Delta_P_inj_old);


d_inj_ox_low = b_d_inj_ox ;
d_inj_f_low = b_d_inj_f;

Delta_P_inj_old = 0.2*Pc;

A_inj_ox = (pi*d_inj_ox^2)/4;
Delta_P_inj_ox_new = (1/(2*rho_ox))*(m_dot_ox/(A_inj_ox*Cd));

A_inj_f = (pi*d_inj_f^2)/4;
Delta_P_inj_f_new = (1/(2*rho_f))*(m_dot_f/(A_inj_f*Cd));

Pc_new_high = Pc_old - (Delta_P_inj_f_new - Delta_P_inj_old) - (Delta_P_inj_f_new - Delta_P_inj_old);
