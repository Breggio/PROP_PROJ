clc
close all
clear all

%% Section 1 - Nominal condition
% Set data
Ru = 8.314472; % [KJ/kmolK]
g0 = 9.81; % [m/s^2]
eps = 80; % ratio A_exit/a_throat, given input
Pc = 20*1e5; % [Pa]
T = 100; % [N]


% Set CEA parameters
OF = 7.2; % from CEA and literature
Cp_1 = 2.5304; % [kJ/KgK] from CEA
T_cc = 2687.2; % [K] from CEA
Mmol = 21.523; % [Kg/Kmol]
Cp = Cp_1*Mmol; % [kJ/kmolK] % mol based 
% C_f2 = 1.8675; % from CEA but it isn't in Vacuum
C_f2 = 1.9646; %[] % C_f in vacuum
Is2 = 3048.3; % m/s 
Pe = 0.01558*1e5; % Pa
Pe2Pc = Pe/Pc;
Is = Is2/g0; % m/s 

% Output
Cv = Cp - Ru; % Mayer's relation
k = Cp/Cv; 
v_exit = sqrt( 1-(Pe2Pc)^((k-1)/k) ) * sqrt( (2*k/(k-1)) * Ru*1e3/Mmol * T_cc); % m/s 

%% Combustion Chamber design
A_t = T/(Pc*C_f2);                                                          %[m^2]
A_t_cm = A_t*1e4;                                                           %[m^2]
D_t_cm = 2*sqrt(A_t_cm/pi);                                                 %[m]
A_e = A_t*eps;                                                              %[m^2]
A_e_cm = A_e*1e4;
D_e_cm = 2*sqrt(A_e_cm/pi);
m_prop = (T-Pe*A_e)/v_exit;
m_ox = (OF/(1+OF)) * m_prop;
m_fu = m_prop - m_ox;
Mach_cc = 0.1; % assumed b/w 0.2 and 0.4
A_cc= (A_t/Mach_cc*((2/(k+1)) * (1+((k-1)/2*Mach_cc^2)))^((k+1)/(2*(k-1)))); %m^2
A_cc_cm = A_cc*1e4; % cm^2
D_cc_cm = 2*sqrt(A_cc_cm/pi); % cm
L_star = 178; % [cm]
V_cc = L_star*A_t_cm; % cm^3 
L_cc = V_cc/A_cc_cm; % cm


%% Conical nozzle
beta = 45*pi/180; % angolo convergente 
alpha = 15*pi/180; % angolo divergente
lambda = (1+ cos(alpha))/2;
T2D = lambda * T; % [N] Thrust considering the divergence loss 
L_DIV = 0.5*(D_e_cm-D_t_cm)/tan(alpha);
L_CON = 0.5*(D_cc_cm-D_t_cm)/tan(alpha);
L_tot = L_CON + L_DIV; % [cm]


lunghezza_camer_di_spinta = L_tot + L_cc; % cm