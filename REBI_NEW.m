close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Tanks

% CHARACTERIZATION OF THE NOMINAL CONDITION: INITAL ONE
% AND SIZING OF THE ENGINE TANKS AND ORIFICES USING THE NOMINAL PARAMETERS

OF = 7.2; % O/F ratio - from previous computations

m_dot = 0.0327; % [kg/s] Propellants mass flow ratio - from previous calculations

m_dot_ox = (OF/(1 + OF))*m_dot; % [kg/s] Oxidizer mass flow rate
m_dot_f = m_dot - m_dot_ox; % [kg/s] Fuel mass flow rate

Pc_in = 20*1e5; % [Pa] Initial combustion chamber pressure

rho_f = 810; % [kg/m^3] fuel density
rho_ox = 1373; % [kg/m^3] oxidizer density

% deltaP_check = 10*6894.76; % [Pa] Pressure loss due to check valve
% deltaP_openclose = 15*6894.76; % [Pa] Pressure loss due to open-close valve
% deltaP_valve = deltaP_check + deltaP_openclose; % Pressure loss due to open-close valve and check valve
% deltaP_feed = 0.05*101325; % [Pa] Pressure loss of the feeding line
d_pipe_f = 0.005; % [m] Diameter of the pipes
A_pipe_f = (d_pipe_f^2)*pi/4; % [m^2] Area of the pipes
d_pipe_ox = 0.005; % [m] Diameter of the oxidizer pipes
A_pipe_ox = (d_pipe_ox^2)*pi/4; % [m^2] Area of the fuel pipes

% [m/s] Velocity of the propellant in the feeding lines
u_f_pipe = m_dot_f/(A_pipe_f*rho_f);
u_ox_pipe = m_dot_ox/(A_pipe_ox*rho_ox);

% deltaP_dyn_f = 0.05*rho_f*u_f^2; % [Pa] Dynamic pressure loss in the feeding lines - fuel
% deltaP_dyn_ox = 0.05*rho_ox*u_ox^2; % [Pa] Dynamic pressure loss in the feeding lines - oxidizer
% deltaP_inj_in = 0.05*Pc_in; % Pressure loss due to injection (15-25% of Pc)
% 
% % Initial pressure in fuel and oxidizer tanks
% Pt_in_f = Pc_in + deltaP_valve + deltaP_feed + deltaP_dyn_f + deltaP_inj_in;
% Pt_in_ox = Pc_in + deltaP_valve + deltaP_feed + deltaP_dyn_ox + deltaP_inj_in;

B = 1.5; % Blow down ratio [3-4]

% % Final pressure in fuel and oxidizer tanks
% Pt_fin_f = Pt_in_f/B;
% Pt_fin_ox = Pt_in_ox/B;

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

%% Injection plate

Cd = 0.65; % [-] Discharge coefficient, depends on geometry & size of plate DA VERIFICARE

N_f = 1; % SCELTO DA NOI
N_ox = 2; % SCELTO DA NOI

Delta_P_inj = 0.2*Pc_in; %IPOTIZZATA COME PRIMA GUESS

A_f = m_dot_f/(Cd*sqrt(2*Delta_P_inj*rho_f)); % [m^2] Fuel total injection area
A_inj_f = A_f/N_f; % [m^2] Area of 1 fuel injector
d_inj_f = sqrt(A_inj_f*4/pi);

A_ox = m_dot_ox/(Cd*sqrt(2*Delta_P_inj*rho_ox)); % [m^2] Oxidizer total injection area
A_inj_ox = A_ox/N_ox; % [m^2] Area of 1 oxidizer injector
d_inj_ox = sqrt(A_inj_ox*4/pi); % [m] Oxidizer injector diameter

u_ox = Cd*sqrt(2*Delta_P_inj/rho_ox); % [m/s] Oxidizer discharge velocity
u_f = Cd*sqrt(2*Delta_P_inj/rho_f); % [m/s] Fuel discharge velocity

gamma_f = 30; % [deg] Oxidizer injector angle, ASSUMED
gamma_ox = asind((m_dot_f/m_dot_ox)*(u_f/u_ox)*sind(gamma_f)); % [deg] Oxidizer injector angle

%% Iterative process
k = 1.66; %He WE DON'T USE THIS, IS IT CORRECT?????
L = 0.50; %[m]  pipe's length MUST BE LOOKED ON SOURCES

% Darcy friction factor from chart (remember to put it in the appendix) for a turbolent 
% flow with a relative pipe roughness of
rough = 0.0025; %[mm]
d_pipe_f_mm = d_pipe_f*10^3; %[mm]
rel_rough_f = rough/d_pipe_f_mm;
d_pipe_ox_mm = d_pipe_ox*10^3; %[mm]
rel_rough_ox = rough/d_pipe_ox_mm;
Re_f = 5.453*10^6;
Re_ox = 6.691*10^5;
f_f = 0.015;
f_ox = 0.014;

[R_cooling_f] = cooling_losses(f_f, d_pipe_f, L, rho_f, A_pipe_f);
[R_cooling_ox] = cooling_losses(f_ox, d_pipe_ox, L, rho_ox, A_pipe_ox);
[R_valves_f] = valves_losses(rho_f,A_pipe_f);
[R_valves_ox] = valves_losses(rho_ox, A_pipe_ox);
[R_inj_f] = inj_loss_reb(rho_f, A_inj_f, Cd, N_f);
[R_inj_ox] = inj_loss_reb(rho_ox, A_inj_ox, Cd, N_ox);
[R_feed_f] = feeding_losses(f_f, rho_f, L, d_pipe_f);
[R_feed_ox] = feeding_losses(f_ox, rho_ox, L, d_pipe_ox);
[R_dyn_f] = dynamic_losses(rho_f, A_f);
[R_dyn_ox] = dynamic_losses(rho_ox, A_ox); 
R_tot_f = R_inj_f + R_feed_f + R_dyn_f + R_valves_f + R_cooling_f;
R_tot_ox = R_inj_ox + R_feed_ox + R_dyn_ox + R_valves_ox + R_cooling_ox;

Pt_in_f = Pc_in + R_tot_f*m_dot_f^2;
Pt_in_ox = Pc_in + R_tot_ox*m_dot_ox^2;

DP_f = R_inj_f*m_dot_f^2;
DP_ox = R_inj_ox*m_dot_ox^2;


OF_vect = OF;
P_tank_f_vect = Pt_in_f;
P_tank_ox_vect = Pt_in_ox;
Pc_vect = Pc_in;
m_dot_f_vect = m_dot_f;
m_dot_ox_vect = m_dot_ox;
m_dot_vect = m_dot;
c_star_old = 1583.7; %DEVONO DARCELO DA CEA
dt = 1; %[cs]
A_t = 2.6772e-5; %[m^2]
m_dot_old = m_dot;
Pc_old = Pc_in;
OF_old = OF;

for i = 1:dt:(tb*10)

    P_tank_f_new = Pt_in_f*(V_gas_in_f / (sum( m_dot_f_vect(1:i))*(dt/10)/rho_f*9.81 + V_gas_in_f ));
    P_tank_f_vect = [P_tank_f_vect P_tank_f_new];
    P_tank_ox_new = Pt_in_ox*(V_gas_in_ox / (sum(m_dot_ox_vect(1:i))*(dt/10)/rho_ox.*9.81 + V_gas_in_ox));
    P_tank_ox_vect = [P_tank_ox_vect P_tank_ox_new];
%     [outputs] = CEA('problem','rocket','frozen','o/f',OF_old,'case','CEAM-rocket1',...
%     'p,Pa',Pc_old,'supsonic(ae/at)',80,'reactants','fuel','RP-1(L)','C',1,...
%     'H',1.95000,'wt%',100,'t(k)',298.0,'oxid','H2O2(L)','wt%',87.5,...
%     't(k)',350,'oxid','H2O(L)','wt%',12.5,'t(k)',350,...
%     'output','thermochemical','end');
%     c_star_new = outputs.output.froz.cstar(1);
    Pc_new = (m_dot_old*c_star_old)/A_t;
    Pc_vect = [Pc_vect Pc_new];
    [R_inj_f] = injection_losses(rho_f, A_inj_f, Cd, N_f);
    [R_inj_ox] = injection_losses(rho_ox, A_inj_ox, Cd, N_ox);
    [R_feed_f] = feeding_losses(f_f, rho_f, L, d_pipe_f);
    [R_feed_ox] = feeding_losses(f_ox, rho_ox, L, d_pipe_ox);
    [R_dyn_f] = dynamic_losses(rho_f, A_f);
    [R_dyn_ox] = dynamic_losses(rho_ox, A_ox);
    R_tot_fuel = R_inj_f + R_feed_f + R_dyn_f;
    R_tot_oxid = R_inj_ox + R_feed_ox + R_dyn_ox;
    [m_dot_f_new] = mass_flow_rate(P_tank_f_new, Pc_new,  R_tot_fuel);
    m_dot_f_vect = [m_dot_f_vect m_dot_f_new];
    [m_dot_ox_new] = mass_flow_rate(P_tank_ox_new, Pc_new, R_tot_ox);
    m_dot_ox_vect = [m_dot_ox_vect m_dot_ox_new];
    m_dot_new = m_dot_ox_new + m_dot_f_new;
    m_dot_vect = [m_dot_vect m_dot_new];
    OF_new = m_dot_ox_new/m_dot_f_new;
    OF_vect = [OF_vect OF_new];

    m_dot_old = m_dot_new;
    OF_old = OF_new;
    Pc_old = Pc_new;
    i


end

% ct_end = outputs.output.froz.cf_vac(3);
% T_end = A_t * Pc_vect(end) * ct_end

%% LOSSES PLOTS

m_dot_f_vect = m_dot_f_vect(1,1:1000);
m_dot_ox_vect = m_dot_ox_vect(1,1:1000);
P_tank_f_vect = P_tank_f_vect(1,1:1000);
P_tank_ox_vect = P_tank_ox_vect(1,1:1000);

subplot(2,3,1)
plot(linspace(1,tb,1000),R_cooling_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,1000),R_cooling_ox.*m_dot_ox_vect.^2)
title('Cooling Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,2)
plot(linspace(1,tb,1000),R_feed_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,1000),R_feed_ox.*m_dot_ox_vect.^2)
title('Feeding Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,3)
plot(linspace(1,tb,1000),R_dyn_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,1000),R_dyn_ox.*m_dot_ox_vect.^2)
title('Dynamic Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,4)
plot(linspace(1,tb,1000),R_inj_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,1000),R_inj_ox.*m_dot_ox_vect.^2)
title('Injection Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')
subplot(2,3,5)
plot(linspace(1,tb,1000),R_valves_f.*m_dot_f_vect.^2)
hold on
plot(linspace(1,tb,1000),R_valves_ox.*m_dot_ox_vect.^2)
title('Velves Losses')
legend('FU','OX')
xlabel('Time [s]')
ylabel('Pressure Losses [Pa]')

% it is clear how the losses which affect the fuel are lower because it's
% mass flow rate is one order of magnitude lower than the one of the
% oxidizer

