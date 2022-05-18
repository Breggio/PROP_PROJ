
% wall 

K_lin = 10; %w/m*k
t = 0.3 ; 

%hot side
T_c = 2687.2; % K
M_c = 0.1; %[/]
a_c = 1106.9; %m/s 5 sonic velocity in CC
v_c = M_c * a_c; %m/S
cp = 2.5304*1e3; % J/Kg*K
T0 = T_c + 0.5*(v_c^2/cp); %K wall - temperature 

% cold side 
T_cooler = 300; %K
K_steel = 16.3; %w/mK
h_cooler = 1500; %W/mK



q = H*(T0 - T_cooler); 

