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

d_pipe = 0.010; % [m] Diameter of the pipes
A_pipe = d_pipe^2*pi/4; % [m^2] Area of the pipes

% [m/s] Velocity of the propellant in the feeding lines
u_f_pipe = m_dot_f/(A_pipe*rho_f);
u_ox_pipe = m_dot_ox/(A_pipe*rho_ox);

B = 2; % Blow down ratio [3-4]

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


%% INJECTION PLATE CORRETTO
% Injection plate

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

%% FRICTION FACTOR CORRETTO
k = 1.66; %He WE DON'T USE THIS, IS IT CORRECT?????
L = 0.20; %[m]  pipe's length MUST BE LOOKED ON SOURCES

% Darcy friction factor from chart (remember to put it in the appendix) for a turbolent 
% flow with a relative pipe roughness of
rough = 0.0025; %[mm]
d_pipe_f_mm = d_pipe*10^3; %[mm]
rel_rough_f = rough/d_pipe_f_mm;
d_pipe_ox_mm = d_pipe*10^3; %[mm]
rel_rough_ox = rough/d_pipe_ox_mm;
Re_f = 5.453*10^6;
Re_ox = 6.691*10^5;
f_f = 0.015;
f_ox = 0.014;

%% LOSSES CORRETTO
% [R_cooling_f] = cooling_losses(f_f, d_pipe_f, L, rho_f, A_pipe_f); 
[R_cooling_ox] = cooling_losses(f_ox, d_pipe, L, rho_ox, A_pipe);
[R_valves_f] = valves_losses(rho_f,A_pipe);
[R_valves_ox] = valves_losses(rho_ox, A_pipe);
[R_inj_f] = inj_loss_reb(rho_f, A_inj_f, Cd, N_f); %DA VERIFICARE CORRETTEZZA EQUAZIONE SBAGLIATA per MAGGI
[R_inj_ox] = inj_loss_reb(rho_ox, A_inj_ox, Cd, N_ox); %DA VERIFICARE CORRETTEZZA EQUAZIONE SBAGLIATA per MAGGI
[R_feed_f] = feeding_losses(f_f, rho_f, L, d_pipe);
[R_feed_ox] = feeding_losses(f_ox, rho_ox, L, d_pipe);
[R_dyn_f] = dynamic_losses(rho_f, A_f);
[R_dyn_ox] = dynamic_losses(rho_ox, A_ox); 
R_tot_f = R_inj_f + R_feed_f + R_dyn_f + R_valves_f;
R_tot_ox = R_inj_ox + R_feed_ox + R_dyn_ox + R_valves_ox + R_cooling_ox;

Pt_in_f = Pc_in + R_tot_f*m_dot_f^2;
Pt_in_ox = Pc_in + R_tot_ox*m_dot_ox^2;

DP_f = R_inj_f*m_dot_f^2; % VERIFICARE CON Delta_P_inj
DP_ox = R_inj_ox*m_dot_ox^2; %VERIFICARE CON Delta_P_inj

%% ITERATIVE PROCESS

dt = 1; %[s]
tb = tb*1; %[s]
c_star = 1583.7; %DEVONO DARCELO DA CEA???????
A_t = 2.5450e-5; %[m^2]

Pc_vect = zeros(tb*dt,1);
m_dot_f_vect = zeros(tb*dt,1);
m_dot_ox_vect = zeros(tb*dt,1);
m_dot_vect = zeros(tb*dt,1);
OF_vect = zeros(tb*dt,1);
P_tank_f_vect = zeros(tb*dt,1);
P_tank_ox_vect = zeros(tb*dt,1);

Pc_vect(1) = Pc_in;
m_dot_f_vect(1) = m_dot_f;
m_dot_ox_vect(1) = m_dot_ox;
P_tank_f_vect(1) = Pt_in_f;
P_tank_ox_vect(1) = Pt_in_ox;
OF_vect(1) = OF;

options = optimset('Display','off');

for i = 2:dt:tb
    
    sum_m_f = dt*sum(m_dot_f_vect);
    sum_m_ox = dt*sum(m_dot_ox_vect);

    fun_f = @(x_f) R_tot_f*x_f^2 + (c_star/A_t)*x_f - Pt_in_f*(V_gas_in_f/(V_gas_in_f + (dt*x_f + sum_m_f)/rho_f));
    m_dot_f_vect(i) = fsolve(fun_f, m_dot_f, options);

    fun_ox = @(x_ox) R_tot_ox*x_ox^2 + (c_star/A_t)*x_ox - Pt_in_ox*(V_gas_in_ox/(V_gas_in_ox + (dt*x_ox + sum_m_ox)/rho_ox));
    m_dot_ox_vect(i) = fsolve(fun_ox, m_dot_ox,options);

    OF_vect(i) = m_dot_ox_vect(i)/m_dot_f_vect(i);

%     [outputs] = CEA('problem','rocket','frozen','o/f',OF_vect(i),'case','CEAM-rocket1',...
%     'p,Pa',Pc_vect(i-1),'supsonic(ae/at)',80,'reactants','fuel','RP-1(L)','C',1,...
%     'H',1.95000,'wt%',100,'t(k)',298.0,'oxid','H2O2(L)','wt%',87.5,...
%     't(k)',350,'oxid','H2O(L)','wt%',12.5,'t(k)',350,...
%     'output','thermochemical','end');
%     c_star = outputs.output.froz.cstar(1);

    m_dot_vect(i) = m_dot_f_vect(i) + m_dot_ox_vect(i);

    Pc_vect(i) = m_dot_vect(i)*c_star/A_t;

    P_tank_f_vect(i) = Pc_vect(i) + R_tot_f*m_dot_f_vect(i)^2;
    P_tank_ox_vect(i) = Pc_vect(i) + R_tot_ox*m_dot_ox_vect(i)^2;
    i

end

% ct_end = outputs.output.froz.cf_vac(3);
% T_end = A_t * Pc_vect(end) * ct_end
% T_1 = A_t * Pc_vect(1) * ct_end;
% T_2 = A_t * Pc_vect(2) * ct_end;
% 
% %% LOSSES PLOTS
% 
% subplot(2,3,1)
% plot(linspace(1,tb,length(m_dot_f_vect)),R_cooling_ox.*m_dot_ox_vect.^2)
% title('Cooling Losses')
% legend('OX')
% xlabel('Time [s]')
% ylabel('Pressure Losses [Pa]')
% subplot(2,3,2)
% plot(linspace(1,tb,length(m_dot_f_vect)),R_feed_f.*m_dot_f_vect.^2)
% hold on
% plot(linspace(1,tb,length(m_dot_f_vect)),R_feed_ox.*m_dot_ox_vect.^2)
% title('Feeding Losses')
% legend('FU','OX')
% xlabel('Time [s]')
% ylabel('Pressure Losses [Pa]')
% subplot(2,3,3)
% plot(linspace(1,tb,length(m_dot_f_vect)),R_dyn_f.*m_dot_f_vect.^2)
% hold on
% plot(linspace(1,tb,length(m_dot_f_vect)),R_dyn_ox.*m_dot_ox_vect.^2)
% title('Dynamic Losses')
% legend('FU','OX')
% xlabel('Time [s]')
% ylabel('Pressure Losses [Pa]')
% subplot(2,3,4)
% plot(linspace(1,tb,length(m_dot_f_vect)),R_inj_f.*m_dot_f_vect.^2)
% hold on
% plot(linspace(1,tb,length(m_dot_f_vect)),R_inj_ox.*m_dot_ox_vect.^2)
% title('Injection Losses')
% legend('FU','OX')
% xlabel('Time [s]')
% ylabel('Pressure Losses [Pa]')
% subplot(2,3,5)
% plot(linspace(1,tb,length(m_dot_f_vect)),R_valves_f.*m_dot_f_vect.^2)
% hold on
% plot(linspace(1,tb,length(m_dot_f_vect)),R_valves_ox.*m_dot_ox_vect.^2)
% title('Velves Losses')
% legend('FU','OX')
% xlabel('Time [s]')
% ylabel('Pressure Losses [Pa]')
% 
% % it is clear how the losses which affect the fuel are lower because it's
% % mass flow rate is one order of magnitude lower than the one of the
% % oxidizer
% 
