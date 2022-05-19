function x = calculate_mdot(x, dt, V_gas_in_ox, V_gas_in_f, C_f, C_ox, At, R_tot_ox, R_tot_f) 

P_tank_ox = Pt_in_ox(V_gas_in_ox/((dt*x(1)+C_ox)+V_gas_in_ox));
P_tank_f = Pt_in_f(V_gas_in_f/((dt*x(2)+C_f)+V_gas_in_f));

[outputs] = CEA('problem','rocket','frozen','o/f',x(1)/x(2),'case','CEAM-rocket1',...
    'p,Pa',x(3),'supsonic(ae/at)',80,'reactants','fuel','RP-1(L)','C',1,...
    'H',1.95000,'wt%',100,'t(k)',298.0,'oxid','H2O2(L)','wt%',87.5,...
    't(k)',350,'oxid','H2O(L)','wt%',12.5,'t(k)',350,...
    'output','thermochemical','end');
c_star = outputs.output.froz.cstar(1);

x(3) = ((x(1)+x(2))*c_star)/At;
X(1) = sqrt((P_tank_ox-x(3))/R_tot_ox);
X(2) = sqrt((P_tank_f-x(3))/R_tot_f);
Pc = ((x(1)+x(2))*c_star)/At;

end
