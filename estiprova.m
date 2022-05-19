





guess_vec = [m_dot_ox, m_dot_f, Pc_in];

fun = @(x) calculate_mdot(dt, V_gas_in_ox, V_gas_in_f, C_f, C_ox, At);

var = fsolve(fun, guess_vec, options);



