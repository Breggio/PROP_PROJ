function F = root3d(x, sum_m_f, sum_m_ox, R_tot_f, R_tot_ox, Pt_in_f, Pt_in_ox, V_gas_in_f, V_gas_in_ox, rho_f, rho_ox, c_star, A_t, dt )
F(1) = R_tot_f*x(1)^2 + x(3) - Pt_in_f*(V_gas_in_f/(V_gas_in_f + (dt/100*x(1) + sum_m_f)/rho_f));
F(2) = R_tot_ox*x(2)^2 + x(3) - Pt_in_ox*(V_gas_in_ox/(V_gas_in_ox + (dt/100*x(2) + sum_m_ox)/rho_ox));
F(3) = x(3) - (x(1) + x(2))*c_star/A_t;
end