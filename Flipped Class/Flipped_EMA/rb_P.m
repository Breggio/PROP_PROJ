function [P_eff_vect, rb_vect] = rb_P(ptrace, web, dt)

%Define geometrical properties
%ptrace [Nx3] [plow, pmid, phigh]
%web   %Propellant thickness, [m]
%dt    %SAmpling time, [s]
% OUTPUT rb mm/s

plow = ptrace(:,1);
pmid = ptrace(:,2);
phigh = ptrace(:,3);

%% Low trace
%Pressure for definition of action time
[Pmax, imax] = max(plow);

if imax  <300 
    imax = round(length(plow)/2);
end

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

I1 = trapz(plow(i_A:i_G)) / 2; %[bar * ms]

%P_ref = I1/(i_G-i_A);
P_ref = I1/((t_G-t_A)*1000);


i = imax;
while (plow(i)>P_ref)
    i = i - 1;
    if i == 0
        break
    end
end
i_B = i;


i = imax;
while (plow(i)>P_ref)
    i = i + 1;

    if i == length(plow)
        break
    end

end
i_E = i;


%t_burn = i_E - i_B;
t_burn_low = (i_E - i_B)*dt; %[s]
P_eff_low = trapz(plow(i_B:i_E)) / (t_burn_low*1000);

rb_low = web/t_burn_low;    %[m/s]

%% Medium trace
%Pressure for definition of action time
[Pmax, imax] = max(pmid);

if imax  <300 
    imax = round(length(pmid)/2);
end

Pact = 0.05 * Pmax;

%To compute action time need to avoid oscillations -> start from Pmax (for
%sure above action time) and move backwards and forward

i = imax;
while (pmid(i)>Pact)
    i = i - 1;
    if i == 0
        break
    end
end
i_A = i;
t_A = dt*(i_A-1);


i = imax;
while (pmid(i)>Pact)
    i = i + 1;

    if i == length(pmid)
        break
    end

end
i_G = i;
t_G = dt*(i_G-1);

I1 = trapz(pmid(i_A:i_G)) / 2; %[bar * ms]

%P_ref = I1/(i_G-i_A);
P_ref = I1/((t_G-t_A)*1000);


i = imax;
while (pmid(i)>P_ref)
    i = i - 1;
    if i == 0
        break
    end
end
i_B = i;


i = imax;
while (pmid(i)>P_ref)
    i = i + 1;

    if i == length(pmid)
        break
    end

end
i_E = i;


%t_burn = i_E - i_B;
t_burn_mid = (i_E - i_B)*dt; %[s]
P_eff_mid = trapz(pmid(i_B:i_E)) / (t_burn_mid*1000);

rb_mid = web/t_burn_mid;    %[m/s]

%% High trace
%Pressure for definition of action time
[Pmax, imax] = max(phigh);

if imax  <300 
    imax = round(length(phigh)/2);
end

Pact = 0.05 * Pmax;

%To compute action time need to avoid oscillations -> start from Pmax (for
%sure above action time) and move backwards and forward

i = imax;
while (phigh(i)>Pact)
    i = i - 1;
    if i == 0
        break
    end
end
i_A = i;
t_A = dt*(i_A-1);


i = imax;
while (phigh(i)>Pact)
    i = i + 1;

    if i == length(phigh)
        break
    end

end
i_G = i;
t_G = dt*(i_G-1);

I1 = trapz(phigh(i_A:i_G)) / 2; %[bar * ms]

%P_ref = I1/(i_G-i_A);
P_ref = I1/((t_G-t_A)*1000);


i = imax;
while (phigh(i)>P_ref)
    i = i - 1;
    if i == 0
        break
    end
end
i_B = i;


i = imax;
while (phigh(i)>P_ref)
    i = i + 1;

    if i == length(phigh)
        break
    end

end
i_E = i;


%t_burn = i_E - i_B;
t_burn_high = (i_E - i_B)*dt; %[s]
P_eff_high = trapz(phigh(i_B:i_E)) / (t_burn_high*1000);

rb_high = web/t_burn_high;    %[m/s]

%%

P_eff_vect = [P_eff_low P_eff_mid P_eff_high];      %[bar]
rb_vect = [rb_low rb_mid rb_high]*1e3;              %[mm/s]

end



