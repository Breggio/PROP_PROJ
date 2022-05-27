function [tb,P] = BARIA(a,n,dt,A_t,c_star)
%INPUT
% a   [mm/s/bar^n]
% A_t [mm^2]
% dt  [s]
% OUTPUT
% P         [bar]
% tb        [s]  

a = a/10^(5*n)/1000; % [m/(s*Pa^n)]
m_ox = 68; %[kg]
m_al = 18;
m_bi = 14;
rho_ox = 1950; % [kg/m^3]
rho_al = 2700;
rho_bi = 920;

rho_av =( m_ox + m_al + m_bi)/(m_ox/rho_ox + m_al/rho_al + m_bi/rho_bi); % [kg/m^3]

l = 290; %  [mm] length of the burning area
x0 = 290/2;
y0 = 30;
dt = 1/1000; %[s]


R_cc = 80; %[mm] radius of c.c



x_new = x0;
y_new = y0;
i = 1;

while x_new >= 0 && y_new >= 0

    Sb = (R_cc-y_new)*2*pi*(2*x_new); %[mm^2]
    Sa = R_cc^2*pi-(R_cc-y_new)^2*pi; %[mm^2]
    A_b = 2*Sa + Sb; %[mm^2]
    P_new = (A_b/A_t*c_star*rho_av*a)^(1/(1-n)); %[Pa]
    P(i) = P_new * 1e-5; %[bar]
    rb = a*P_new^n * 1e3; %[mm/s]
    
    x_old = x_new;
    y_old = y_new;
    x_new = x_old - dt*rb; %[mm]
    y_new = y_old - dt*rb; %[mm]

    i = i + 1;

end

tb = dt*(i-1); %[s]

end