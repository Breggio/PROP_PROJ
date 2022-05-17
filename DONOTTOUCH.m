close all; clearvars; clc;
set(0,'defaulttextInterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Tanks

% CHARACTERIZATION OF THE NOMINAL CONDITION: INITAL ONE
% AND SIZING OF THE ENGINE TANKS AND ORIFICES USING THE NOMINAL PARAMETERS

O_F = 6; % O/F ratio - from previous computations

m_dot = 0.0312; % [kg/s] Propellants mass flow ratio - from previous calculations

m_dot_ox = (O_F/(1 + O_F))*m_dot; % [kg/s] Oxidizer mass flow rate
m_dot_f = m_dot - m_dot_ox; % [kg/s] Fuel mass flow rate

Pc_in = 20*1e5; % [Pa] Initial combustion chamber pressure

rho_f = 810; % [kg/m^3] fuel density
rho_ox = 1373; % [kg/m^3] oxidizer density

deltaP_check = 10*6894.76; % [Pa] Pressure loss due to check valve
deltaP_openclose = 15*6894.76; % [Pa] Pressure loss due to open-close valve
deltaP_valve = deltaP_check + deltaP_openclose; % Pressure loss due to open-close valve and check valve
deltaP_feed = 0.5*101325; % [Pa] Pressure loss of the feeding line
d_pipe = 0.005; % [m] Diameter of the pipes
A_pipe = d_pipe^2*pi/4; % [m^2] Area of the pipes

% [m/s] Velocity of the propellant in the feeding lines
u_f = m_dot_f/A_pipe/rho_f;
u_ox = m_dot_ox/A_pipe/rho_ox;

deltaP_dyn_f = 0.5*rho_f*u_f^2; % [Pa] Dynamic pressure loss in the feeding lines - fuel
deltaP_dyn_ox = 0.5*rho_ox*u_ox^2; % [Pa] Dynamic pressure loss in the feeding lines - oxidizer
deltaP_inj_in = 0.2*Pc_in; % Pressure loss due to injection (15-25% of Pc)

% Initial pressure in fuel and oxidizer tanks
Pt_in_f = Pc_in + deltaP_valve + deltaP_feed + deltaP_dyn_f + deltaP_inj_in;
Pt_in_ox = Pc_in + deltaP_valve + deltaP_feed + deltaP_dyn_ox + deltaP_inj_in;

B = 3; % Blow down ratio [3-4]

% % Final pressure in fuel and oxidizer tanks
% Pt_fin_f = Pt_in_f/B;
% Pt_fin_ox = Pt_in_ox/B;

tb = 100; % [s] Burning time
M_f = m_dot_f*tb; % [kg] Fuel mass
M_ox = m_dot_ox*tb; % [kg] Oxidizer mass

V_f = M_f*rho_f; % [m^3] Fuel volume
V_ox = M_ox*rho_ox; % [m^3] Oxidizer volume

% Volume occupied initially by the pressurized gas (HYPOTESIS)
V_gas_in_f = V_f/(B-1); % [m^3] Fuel tank
V_gas_in_ox = V_ox/(B-1); % [m^3] Oxidizer tank

% Volume occupied finally by the pressurized gas = tank volume
V_tank_f = V_f + V_gas_in_f; % [m^3] Fuel tank
V_tank_ox = V_ox + V_gas_in_ox; % [m^3] Oxidizer tank

%% Injection plate3

% P_t_f_iso = Pt_in_f.*( V_gas_in_f ./ (m_dot_f.*time./rho_f.*9.81 + V_gas_in_f));


















