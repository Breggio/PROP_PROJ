%% Design Combustion Chamber and Nozzle (Conical)
% Assumption: - Frozen condition across the whole system
% - Storage Temperature Fuel: 298.15K 
% - Storage Temperature Oxidizer: 350K shall be < 423K (boiling)
clc
clear all
close all
%% Section 1 - Nominal condition
% Set data
R = 8.314472; % [KJ/kmolK] or [J/molK]
g0 = 9.807; % [m/s^2]
eps = 80; % ratio A_exit/a_throat, given input
Pc = 20*1e5; % [Pascal], given input
T = 100; % [N], given NOMINAL input
beta = 45*pi/180;
alpha = 15*pi/180;
lambda = (1+cos(alpha))/2; % coff of losses 
T2D = lambda * T; 
T_new = T + T*(1-lambda)
% Set CEA parameters
OF = 7.2; % from CEA and literature
Cp_1 = 2.5322; % [kJ/KgK] from CEA
T_cc = 2699.21; % [K] from CEA
Mmol = 21.506; % [kg/kmol]
Cp = Cp_1*Mmol; % [kJ/kmolK]
C_f2 = 1.8676; % [-]
Is2 = 3056.5; % [m/s]
Pe = 0.01559*10^5; % [Pascal]
c_star = 1583.7; % [m/s]
Pe2Pc = Pe/Pc;
Is = Is2/g0; % [s]
% Output
Cv = Cp - R; % Mayer's relation
k = Cp/Cv; % definition of gamma/k, and also from CEA
Mach_e = 4.856; % from CEA
v_sonic_e = 609.1; % [m/s]
u_e = Mach_e * v_sonic_e; % exit velocity [m/s]

% Chamber design
A_t = T_new/(Pc*C_f2); % m 
A_t_cm = A_t*1e4;
D_t_cm = 2*sqrt(A_t_cm/pi);
A_e = A_t*eps; % m
% m_dot2 = Pc*A_t/c_star;
m_dot = (T_new-Pe*A_e)/u_e;
m_ox = (OF/(1+OF)) * m_dot;
m_fu = m_dot - m_ox;
A_e_cm = A_e*10^4;
D_e_cm = 2*sqrt(A_e_cm/pi);
%m_dot = (T-Pe*A_e)/v_exit;
Mach_cc = 0.1; % assumed b/w 0.2 and 0.4
A_cc_cm = (A_t/Mach_cc*((2/(k+1)) * (1+((k-1)/2*Mach_cc^2)))^((k+1)/2/(k-1)))*1e4;
D_cc_cm = 2*sqrt(A_cc_cm/pi);
L_star = 178*1.1; % [cm]
V_cc = L_star*A_t_cm; % cm^3
L_cc = (V_cc/A_cc_cm); % cm
contract_ratio = A_cc_cm / A_t_cm
%% Conical nozzle
% beta = 45*pi/180;
% alpha = 15*pi/180;
% lambda = (1+cos(alpha))/2; % coff of losses 
% T2D = lambda * T_new; 
L_DIV = 0.5*(D_e_cm-D_t_cm)/tan(alpha); % cm
L_CON = 0.5*(D_cc_cm-D_t_cm)/tan(beta); % cm
L_tot_nozzle = L_CON + L_DIV; % cm

L_engine = L_cc + L_tot_nozzle
