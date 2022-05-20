%%

close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Design Combustion Chamber and Nozzle (Conical)

% Assumption: - Frozen condition across the whole system
% - Storage Temperature Fuel: 298.15K 
% - Storage Temperature Oxidizer: 350K shall be < 423K (boiling)

%% Section 1 - Nominal condition

% Set data

R = 8.314472; % [KJ/kmolK] or [J/molK]
g0 = 9.807; % [m/s^2]
eps = 80; % ratio A_exit/a_throat, given input
Pc = 20*1e5; % [Pascal], given input
T = 100; % [N], GIVEN NOMINAL input
beta = 45*pi/180; % SUPPOSED
alpha = 15*pi/180; % SUPPOSED
lambda = (1+cos(alpha))/2; % coff of losses 
T = (2-lambda) * T;  %NEW NOMINAL THURST WHICH WILL BE AFFECTED BY THE LOSSES

% Set CEA parameters
OF = 7.2; % from CEA and literature
Cp_1 = 2.5322; % [kJ/KgK] from CEA
T_cc = 2699.21; % [K] from CEA
Mmol = 21.506; % [kg/kmol]
Cp = Cp_1*Mmol; % [kJ/kmolK]
% C_f = 1.8676; % [-]
Is2 = 3056.5; % [m/s]
Pe = 0.01559*10^5; % [Pascal]
% c_star = 1583.7; % [m/s]
Pe2Pc = Pe/Pc;
Is = Is2/g0; % [s]
% Output
Cv = Cp - R; % Mayer's relation [kJ/kmolK]
k = Cp/Cv; % definition of gamma/k, and also from CEA
Mach_e = 4.856; % from CEA
v_sonic_e = 609.1; % [m/s]
u_e = Mach_e * v_sonic_e; % exit velocity [m/s]
C_f = sqrt(2*(k^2/(k-1))*(2/(k+1))^((k+1)/(k-1)))*sqrt(1-(Pe2Pc)^((k-1)/k))+eps*Pe2Pc;
c_star = sqrt(k*(R*1e3/Mmol)*T_cc)*1/(k*sqrt((2/(k+1))^((k+1)/(k-1))));

% Chamber design
A_t = T/(Pc*C_f); % m 
A_t_cm = A_t*1e4;
D_t_cm = 2*sqrt(A_t_cm/pi);
A_e = A_t*eps; % m
m_dot = Pc*A_t/c_star;
m_dot_ox = (OF/(1+OF)) * m_dot;
m_dot_f = m_dot - m_dot_ox;
A_e_cm = A_e*10^4;
D_e_cm = 2*sqrt(A_e_cm/pi);
Mach_cc = 0.1; % assumed below 0.6
A_cc_cm = (A_t/Mach_cc*((2/(k+1)) * (1+((k-1)/2*Mach_cc^2)))^((k+1)/2/(k-1)))*1e4;
D_cc_cm = 2*sqrt(A_cc_cm/pi);
L_star = 178*1.1; % [cm]
V_cc = L_star*A_t_cm; % cm^3
L_cc = (V_cc/A_cc_cm); % cm
contract_ratio = A_cc_cm/A_t_cm;

%% Conical nozzle

L_DIV = 0.5*(D_e_cm-D_t_cm)/tan(alpha); % cm
L_CON = 0.5*(D_cc_cm-D_t_cm)/tan(beta); % cm
L_tot_nozzle = L_CON + L_DIV; % cm

L_engine = L_cc + L_tot_nozzle;


%% BLOW-DOWN DESIGN
%% Tanks

% CHARACTERIZATION OF THE NOMINAL CONDITION: INITAL ONE
% AND SIZING OF THE ENGINE TANKS AND ORIFICES USING THE NOMINAL PARAMETERS

rho_f = 810; % [kg/m^3] fuel density
rho_ox = 1373; % [kg/m^3] oxidizer density

d_pipe_f = 0.005; % [m] Diameter of the pipes
A_pipe_f = d_pipe_f^2*pi/4; % [m^2] Area of the pipes
d_pipe_ox = 0.010; % [m] Diameter of the pipes
A_pipe_ox = d_pipe_ox^2*pi/4; % [m^2] Area of the pipes

% [m/s] Velocity of the propellant in the feeding lines
u_f_pipe = m_dot_f/(A_pipe_f*rho_f);
u_ox_pipe = m_dot_ox/(A_pipe_ox*rho_ox);

B = 1.2; % Blow down ratio [3-4]

% % Final pressure in fuel and oxidizer tanks
tb = 100; % [s] Burning time

M_f = m_dot_f*tb; % [kg] Fuel mass
M_ox = m_dot_ox*tb; % [kg] Oxidizer mass

V_f = M_f/rho_f; % [m^3] Fuel volume
V_ox = M_ox/rho_ox; % [m^3] Oxidizer volume

% Volume occupied initially by the pressurized gas (HYPOTESIS)
V_gas_in_f = V_f/(B-1); % [m^3] Fuel tank
V_gas_in_ox = V_ox/(B-1); % [m^3] Oxidizer tank

% Volume occupied finally by the pressurized gas = tank volume
V_tank_f = V_f + V_gas_in_f; % [m^3] Fuel tank
V_tank_ox = V_ox + V_gas_in_ox; % [m^3] Oxidizer tank
V_tank_f = V_tank_f * 1.05; % [m^3] Fuel tank + margins 
V_tank_ox = V_tank_ox * 1.05; % [m^3] Oxidizer tank + margins

%% INJECTION PLATE CORRETTO
% Injection plate

Cd = 0.7; % [-] Discharge coefficient, depends on geometry & size of plate DA VERIFICARE

N_f = 1; % SCELTO DA NOI
N_ox = 2; % SCELTO DA NOI

Delta_P_inj = 0.20*Pc; %IPOTIZZATA COME PRIMA GUESS

A_f = m_dot_f/(Cd*sqrt(2*Delta_P_inj*rho_f)); % [m^2] Fuel total injection area
A_inj_f = A_f/N_f; % [m^2] Area of 1 fuel injector
d_inj_f = sqrt(A_inj_f*4/pi);

A_ox = m_dot_ox/(Cd*sqrt(2*Delta_P_inj*rho_ox)); % [m^2] Oxidizer total injection area
A_inj_ox = A_ox/N_ox; % [m^2] Area of 1 oxidizer injector
d_inj_ox = sqrt(A_inj_ox*4/pi); % [m] Oxidizer injector diameter

u_ox = Cd*sqrt(2*Delta_P_inj/rho_ox); % [m/s] Oxidizer discharge velocity
u_f = Cd*sqrt(2*Delta_P_inj/rho_f); % [m/s] Fuel discharge velocity

%% FRICTION FACTOR CORRETTO

L_pipes = 0.50; %[m]  pipe's length MUST BE LOOKED ON SOURCES
n_cool = 15;
D_cc_m = 1.4261*1e-2;
L_cool = n_cool*D_cc_m*pi;

% Darcy friction factor from chart (remember to put it in the appendix) for a turbolent 
% flow with a relative pipe roughness of

rough = 0.0025; %[mm]
d_pipe_f_mm = d_pipe_f*10^3; %[mm]
rel_rough_f = rough/d_pipe_f_mm;
d_pipe_ox_mm = d_pipe_ox*10^3; %[mm]
rel_rough_ox = rough/d_pipe_ox_mm;
mu_ox = 1.26e-3; %pa*s 
mu_f = 1.40e-3; %pa*s
Re_ox = (rho_ox * d_pipe_ox * u_ox_pipe)/mu_ox;
Re_f = (rho_f * d_pipe_f * u_f_pipe)/mu_f; 
f_f = 64/Re_f; % FOR LAMINAR FLOW
f_ox = 0.045; % FROM CHART

%% LOSSES CORRETTO

[R_cooling_ox] = cooling_losses(f_ox, d_pipe_ox, L_cool, rho_ox, A_pipe_ox);
[R_valves_f] = control_valve_pressure_loss(rho_f);
[R_valves_ox] = control_valve_pressure_loss(rho_ox);
[R_inj_f] = injection_losses(rho_f, A_inj_f, Cd, N_f); 
[R_inj_ox] = injection_losses(rho_ox, A_inj_ox, Cd, N_ox); 
[R_feed_f] = feeding_losses(f_f, rho_f, L_pipes, d_pipe_f);
[R_feed_ox] = feeding_losses(f_ox, rho_ox, L_pipes, d_pipe_ox);
[R_dyn_f] = dynamic_losses(rho_f, A_f);
[R_dyn_ox] = dynamic_losses(rho_ox, A_ox); 
R_tot_f = R_inj_f + R_feed_f + R_dyn_f + R_valves_f;
R_tot_ox = R_inj_ox + R_feed_ox + R_dyn_ox + R_valves_ox + R_cooling_ox;

Pt_in_f = Pc + R_tot_f*m_dot_f^2;
Pt_in_ox = Pc + R_tot_ox*m_dot_ox^2;

DP_f = R_inj_f*m_dot_f^2; 
DP_ox = R_inj_ox*m_dot_ox^2; 

%% ITERATIVE PROCESS

dt = 1; %[s]
tb = tb*1000; %[s]

Pc_vect = zeros(1,tb*dt);
m_dot_f_vect = zeros(1,tb*dt);
m_dot_ox_vect = zeros(1,tb*dt);
m_dot_vect = zeros(1,tb*dt);
OF_vect = zeros(1,tb*dt);
sum_m_f = zeros(1,tb*dt);
sum_m_ox = zeros(1,tb*dt);
T_vect = zeros(1,tb*dt);

Pc_vect(1) = Pc;
m_dot_f_vect(1) = m_dot_f;
m_dot_ox_vect(1) = m_dot_ox;
OF_vect(1) = OF;
T_vect(1) = T;

options = optimset('Display','off');

for i = 2:dt:tb

    sum_m_f(i-1) = dt*sum(m_dot_f_vect);
    sum_m_ox(i-1) = dt*sum(m_dot_ox_vect);

    fun = @(x)root3d(x, sum_m_f(i-1), sum_m_ox(i-1), R_tot_f, R_tot_ox, Pt_in_f, Pt_in_ox, V_gas_in_f, V_gas_in_ox, rho_f, rho_ox, c_star, A_t, dt);
    x0 = [m_dot_f_vect(i-1), m_dot_ox_vect(i-1), Pc_vect(i-1)];
    x = fsolve(fun, x0, options);

    m_dot_f_vect(i) = x(1);
    m_dot_ox_vect(i) = x(2);
    Pc_vect(i) = x(3);

    OF_vect(i) = m_dot_ox_vect(i)/m_dot_f_vect(i);

    p_c_star = CEA_int_C_star();
    c_star = Val_Int_C_star(p_c_star,Pc_vect(i)*1e-5,OF_vect(i));
    p_ct = CEA_int_C_F(A_e, A_t);
    ct = Val_Int_C_F(p_ct,Pc_vect(i)*1e-5,OF_vect(i));

    T_vect(i) = A_t * Pc_vect(i) * ct;

end

P_tank_f_vect = Pc_vect + R_tot_f.*m_dot_f_vect.^2;
P_tank_ox_vect = Pc_vect + R_tot_ox.*m_dot_ox_vect.^2;


%% LOSSES PLOTS

subplot(2,3,1)
plot(linspace(1,tb,length(m_dot_f_vect)),R_cooling_ox.*m_dot_ox_vect.^2)
title('Cooling Losses')
legend('OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,2)
plot(linspace(1,tb,length(m_dot_f_vect)),R_feed_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,length(m_dot_f_vect)),R_feed_ox.*m_dot_ox_vect.^2)
title('Feeding Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,3)
plot(linspace(1,tb,length(m_dot_f_vect)),R_dyn_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,length(m_dot_f_vect)),R_dyn_ox.*m_dot_ox_vect.^2)
title('Dynamic Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,4)
plot(linspace(1,tb,length(m_dot_f_vect)),R_inj_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,length(m_dot_f_vect)),R_inj_ox.*m_dot_ox_vect.^2)
title('Injection Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,5)
plot(linspace(1,tb,length(m_dot_f_vect)),R_valves_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,length(m_dot_f_vect)),R_valves_ox.*m_dot_ox_vect.^2)
title('Velves Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')

% it is clear how the losses which affect the fuel are lower because it's
% mass flow rate is one order of magnitude lower than the one of the
% oxidizer

%% PLOT

time = [1:1:tb];

%combustion chamber pressure
figure()
plot(time, Pc_vect,'c', 'LineWidth', 2.5)
legend('Combustion chamber pressure', 'FontSize', 30)
grid on; grid minor
xlabel('Time [s]', 'FontSize', 30)
ylabel('Pc [Pa]', 'FontSize', 30)

%thrust
figure()
plot(time, T_vect,'m', 'LineWidth', 2.5)
legend('Thrust', 'FontSize', 30)
grid on; grid minor
xlabel('Time [s]', 'FontSize', 30)
ylabel('T [N]', 'FontSize', 30)

%specific impulse
g_0 = 9.81;
m_dot_vect = m_dot_ox_vect + m_dot_f_vect;
I_sp_vec = T_vect./ (m_dot_vect.* g_0);
figure()
plot(time, I_sp_vec, 'g', 'LineWidth', 2.5)
legend('Specific impulse', 'FontSize', 30)
grid on; grid minor
xlabel('Time [s]', 'FontSize', 30)
ylabel('$I_{sp} [s]$','FontSize', 30, 'Interpreter', 'latex')

%m_dot
% figure()
% plot(time, m_dot_ox_vect, 'g', 'LineWidth', 2.5)
% hold on 
% plot(time,  m_dot_f_vect, 'k', 'LineWidth', 2.5)
% legend('$\dot{m_{ox}}$', '$\dot{m_{f}}$', 'FontSize', 30)
% grid on; grid minor
% xlabel('$Time [s]$', 'FontSize', 30, 'Interpreter', 'latex')
% ylabel('$\dot{m} [kg]$','FontSize', 30, 'Interpreter', 'latex')

figure()
plot(time, m_dot_ox_vect, 'g', 'LineWidth', 2.5)
% legend('$\dot{m_{ox}}$', 'FontSize', 30)
grid on; grid minor
xlabel('$Time [s]$', 'FontSize', 30, 'Interpreter', 'latex')
ylabel('$\dot{m_{ox}} [kg]$','FontSize', 30, 'Interpreter', 'latex')

%m_fuel
figure()
plot(time,  m_dot_f_vect, 'g', 'LineWidth', 2.5)
% legend('Fuel mass flow rate', 'FontSize', 30)
grid on; grid minor
xlabel('$Time [s]$', 'FontSize', 30)
ylabel('$\dot{m_{f}} [kg]$','FontSize', 30, 'Interpreter', 'latex')

%OF
figure()
plot(time,  OF_vect, 'b', 'LineWidth', 2.5)
legend('$OF$', 'FontSize', 30, 'Interpreter', 'latex')
grid on; grid minor
xlabel('$Time [s]$', 'FontSize', 30)
ylabel('$OF [-]$','FontSize', 30, 'Interpreter', 'latex')

%Pressure tank ox and f
% figure()
% plot(time, P_tank_ox_vect, 'k', 'LineWidth', 2.5)
% hold on
% plot(time, P_tank_f_vect, 'c', 'LineWidth', 2.5)
% legend('$P_{t,ox}$','$P_{t,f}$', 'FontSize', 30, 'Interpreter', 'latex')
% grid on; grid minor
% xlabel('$Time [s]$', 'FontSize', 30)
% ylabel('$Pressure [Pa]$','FontSize', 30, 'Interpreter', 'latex')
