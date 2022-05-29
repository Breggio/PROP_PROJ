function [tb, P] = BARIA(a, n, delta_t, At, c_star, rho_p, Time)

x_old = 0.29/2; % [m]
y_old = 0.03; % [m]

S_A = (0.08)^2*pi - (0.08-y_old)^2*pi;
S_B = (0.08-y_old)*2*pi*2*x_old;

Ab = 2*S_A + S_B; % [m^2]
a = a*10^(-5*n);
P_in = ((Ab/At)*c_star*rho_p*a)^(1/(1-n)); % [bar]
Rb = a*(P_in)^n;

for i=1:Time

x_new = x_old - delta_t*Rb;
y_new = y_old - delta_t*Rb;

S_A = (0.08)^2*pi - (0.08-y_new)^2*pi;
S_B = (0.08-y_new)*2*pi*2*x_new;

Ab = 2*S_A + S_B;

P(i) = ((Ab/At)*c_star*rho_p*a)^(1/(1-n));
Rb = a*(P(i))^n;

x_old = x_new; 
y_old = y_new; 

if x_new <= 0 || y_new <= 0
    tb = i/1e3;
    break
end

end

end