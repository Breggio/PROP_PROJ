function [tb, P_B] = BARIA(a_Pa, n, delta_t, At, c_star, rho_p, Time)

% Inputs 
% a_Pa [mm/(s Pa^n)]
% ...
%
% Outputs
%
% tb [ms]
% P [Pa]

x_old = 290/2; % [mm]
y_old = 30; % [mm]
r_cc = 80; % [mm] c.c radius

S_A = r_cc^2*pi - (r_cc-y_old)^2*pi; % [mm^2]
S_B = (r_cc-y_old)*2*pi*2*x_old; % [mm^2]

Ab = 2*S_A + S_B; % [mm^2]

P_in = ((Ab/At)*c_star*rho_p*a_Pa)^(1/(1-n)); % [Pa]
Rb = a_Pa*(P_in/1e5)^n;


for i=1:Time

    x_new = x_old - delta_t*Rb;
    y_new = y_old - delta_t*Rb;

    S_A = (r_cc)^2*pi - (r_cc-y_new)^2*pi;
    S_B = (r_cc-y_new)*2*pi*2*x_new;

    Ab = 2*S_A + S_B;
    
    P(i) = ( ((Ab/At)*c_star*rho_p*a_Pa)^(1/(1-n)) ); % [Pa]
    Rb = a_Pa*(P(i))^n;
    
    x_old = x_new; 
    y_old = y_new; 
    
    if x_new <= 0 || y_new <= 0
    tb = i; % [s]
    disp('COGLIONE')
    break
    end

end
P_B = P*1e-5;
end