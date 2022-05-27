clc
clear all
close all

%Define geometrical properties
web = 0.03;     %Propellant thickness, 3 cm
dt = 1/1000;    %SAmpling time, [s]

load("tracesbar1.mat");

plow = pbar2438(:,1);
pmid = pbar2438(:,2);
phigh = pbar2438(:,3);

%Pressure for definition of action time
[Pmax, imax] = max(plow);
Pact = 0.05 * Pmax;

%To compute action time need to avoid oscillations -> start from Pmax (for
%sure above action time) and move backwards and forward

i = imax;
while (plow(i)>Pact)
    i = i - 1;
    if i == 0
        break
    end
end
i_A = i;
t_A = dt*(i_A-1);


i = imax;
while (plow(i)>Pact)
    i = i + 1;

    if i == length(plow)
        break
    end

end

i_G = i;
t_G = dt*(i_G-1);

I1 = trapz(plow(i_A:i_G)) / 2;

P_ref = I1/(i_G-i_A)
%P_ref = I1/(t_G-t_A)
%%

i = imax;

while (plow(i)>P_ref)
    i = i - 1;
    if i == 0
        break
    end
end
i_B = i;
t_B = dt*(i_B-1);

i = imax;

while (plow(i)>P_ref)
    i = i + 1;

    if i == length(plow)
        break
    end

end

i_E = i;


%t_burn = i_E - i_B;
t_burn = (i_E - i_B)*dt;
P_eff = trapz(plow(i_B:i_E)) / t_burn

rb = web/t_burn;    %[m/s]
