clear
clc

load('tracesbar1.mat')

%Define geometrical properties
web = 0.03;     %Propellant thickness [m]
dt = 1/1000;    %SAmpling time, [s]

% CHECK ERROR IN DIMENSION NON HO CAPITO CHE CAZZO VOLESSE DIRE, CHE ERRORE
% C'E??

[P_eff_2438, rb_2438] = rb_P(pbar2438, web, dt);
[P_eff_2439, rb_2439] = rb_P(pbar2439, web, dt);
[P_eff_2440, rb_2440] = rb_P(pbar2440, web, dt);
[P_eff_2441, rb_2441] = rb_P(pbar2441, web, dt);
[P_eff_2442, rb_2442] = rb_P(pbar2442, web, dt);
[P_eff_2443, rb_2443] = rb_P(pbar2443, web, dt);
[P_eff_2444, rb_2444] = rb_P(pbar2444, web, dt); 
%In 2444 high pressure inverted with inetermediate pressure trace ->
%reorder them
P_eff_2444 = [P_eff_2444(1),P_eff_2444(3),P_eff_2444(2)];
rb_2444 = [rb_2444(1),rb_2444(3),rb_2444(2)];
[P_eff_2445, rb_2445] = rb_P(pbar2445, web, dt);
[P_eff_2446, rb_2446] = rb_P(pbar2446, web, dt);


P_eff_vect = [P_eff_2438,P_eff_2439,P_eff_2440,P_eff_2441,P_eff_2442,P_eff_2443,P_eff_2444,P_eff_2445,P_eff_2446];%[bar]
rb_vect = [rb_2438,rb_2439,rb_2440,rb_2441,rb_2442,rb_2443,rb_2444,rb_2445,rb_2446]; %[mm/s]
%%
figure
hold on
plot(pbar2438(:,1),'LineWidth',1.5)
plot(pbar2438(:,2),'LineWidth',1.5)
plot(pbar2438(:,3),'LineWidth',1.5)
legend('Low pressure','Mid pressure','High pressure')
title('Pressure traces')
xlabel('Time [s]')
ylabel('Pressure [bar]')
grid on
%%
% Plot
figure
loglog(P_eff_vect,rb_vect, 'o','LineWidth',1);
grid on
hold on

% Calcolo parametri vieille law
[a, Inc_a, n, Inc_n, R2] = Uncertainty(P_eff_vect, rb_vect);
% a [mm/s/bar^n]
% Inc_a [mm/s/bar^n]

rb_fun = @(P) a.*P.^n;
loglog(linspace(P_eff_2446(1)-1,P_eff_2446(3)+1,100), rb_fun(linspace(P_eff_2446(1)-1,P_eff_2446(3)+1,100)),'LineWidth',1.5)
title('Fitting of the Vieille''s law' )
legend('Experimental data','Fitting Vieille''s law')
xlabel('Pressure [bar]')
ylabel('Burning rate [mm/s]')

%% MONTECARLO

N = 100;



a_pop = normrnd(a,Inc_a,[1,N]);
n_pop = normrnd(n,Inc_n,[1,N]);
tb_vect = [];
tb_mean = [];
tb_dev = [];

[A_pop, N_pop] = meshgrid(a_pop,n_pop);

[x,y] = size(A_pop);

%%
for i = 1:x
    %randperm(x) crea un vettore di numeri da 1 a x ordinati random
    idx1 = randperm(x);
    A_pop(i,:) = A_pop(i,idx1);
    %la riga i-esima è uguale a quella di prima con le colonne riordinate
    %randomicamente
    N_pop(i,:) = N_pop(i,idx1); %NON NECESSARIO SICCOME NUmERI UGUALI PER RIGA -> COPPIE SONO MANTENUTE
    %In questo modo ho rioridnato randomicamente le coppie movendole sulla
    %riga
   
end

 %devo fare shuffle anche su colonna della stessa matrice
for j = 1:y

        idx2 = randperm(y);
        A_pop(:,j) = A_pop(idx2,j);
        %la riga i-esima è uguale a quella di prima con le colonne riordinate
        %randomicamente
        N_pop(:,j) = N_pop(idx2,j);
        %In questo modo ho rioridnato randomicamente le coppie movendole sulla
        %colonna
end

%%
% 
% for i = 1:x
%     
% %     idx1 = randperm(x);
% %     A_pop(i,:) = A_pop(i,idx1);
%     %N_pop(i,:) = N_pop(i,idx1); %NON NECESSARIO SICCOME NUmERI UGUALI PER RIGA -> COPPIE SONO MANTENUTE
%     %In questo modo ho rioridnato randomicamente le coppie movendole sulla
%     %riga
% 
%     for j = 1:y
% 
%         idx2 = randperm(y);
%         A_pop(:,j) = A_pop(idx2,j);
%         N_pop(:,j) = N_pop(idx2,j);
% 
%         for k = 1:x
%             idx3 = randperm(x);
%             A_pop(k,:) = A_pop(k,idx3);
%             N_pop(k,:) = N_pop(k,idx3);
%         end
%     end
%    
% end



%% LOW
tb_vect_low = [];
tb_mean_low = [];
tb_dev_low = [];
err_t = [];
err_dev = [];

A_t_low = pi*(28.8/2)^2; %[mm^2]
c_star_low = 1610.4;

for i = [1:N] 

    for j = [1:N]

        %[tb,P] = BARIA(a_pop(i),n_pop(j),dt,A_t_low); %no shuffle
        [tb,P] = BARIA(A_pop(i,j),N_pop(i,j),dt,A_t_low,c_star_low); %shuffle righe colonne 2
        tb_vect_low = [tb_vect_low tb];
        tb_mean_low = [tb_mean_low mean(tb_vect_low)];
        tb_dev_low = [tb_dev_low std(tb_vect_low)];
        
        if length(tb_mean_low)>1
            delta_t  = abs(tb_mean_low(end) - tb_mean_low(end-1))/abs(tb_mean_low(end));
            delta_dev = abs(tb_dev_low(end) - tb_dev_low(end-1))/abs(tb_dev_low(end));
            err_t = [err_t delta_t];
            err_dev = [err_dev delta_dev];
        end
    end
end
%%
figure 
semilogy(err_t)
legend('\delta_1')
xlabel('iteration')
figure 
semilogy(err_dev)
legend('\delta_2')
xlabel('iteration')

%%
figure
plot(tb_mean_low)
title('t_B mean')
figure
plot(tb_dev_low)
title('t_B std dev')
tb_low = tb_mean_low(end);
tb_low_std_dev = tb_dev_low(end);

err_tb_mean_low = abs((tb_mean_low(end)-tb_mean_low(end-1))/tb_mean_low(end));
err_tb_dev_low = abs((tb_dev_low(end)-tb_dev_low(end-1))/tb_dev_low(end));

%% MID
tb_vect_mid = [];
tb_mean_mid = [];
tb_dev_mid = [];

A_t_mid = pi*(25.25/2)^2; %[mm^2]
c_star_mid = 1616.6;


for i = [1:N] 

    for j = [1:N]

        %[tb,P] = BARIA(a_pop(i),n_pop(j),dt,A_t_low); %no shuffle
        [tb,P] = BARIA(A_pop(i,j),N_pop(i,j),dt,A_t_mid,c_star_mid); %shuffle righe colonne 2
        tb_vect_mid = [tb_vect_mid tb];
        tb_mean_mid = [tb_mean_mid mean(tb_vect_mid)];
        tb_dev_mid = [tb_dev_mid std(tb_vect_mid)];

    end
end
%%
err_tb_mean_mid = abs((tb_mean_mid(end)-tb_mean_mid(end-1))/tb_mean_mid(end));
err_tb_dev_mid = abs((tb_dev_mid(end)-tb_dev_mid(end-1))/tb_dev_mid(end));
%% HIGH
tb_vect_high = [];
tb_mean_high = [];
tb_dev_high = [];

A_t_high = pi*(21.81/2)^2; %[mm^2]
c_star_high = 1622.8;


for i = [1:N] 

    for j = [1:N]

        %[tb,P] = BARIA(a_pop(i),n_pop(j),dt,A_t_low); %no shuffle
        [tb,P] = BARIA(A_pop(i,j),N_pop(i,j),dt,A_t_high,c_star_high); %shuffle righe colonne 2
        tb_vect_high = [tb_vect_high tb];
        tb_mean_high = [tb_mean_high mean(tb_vect_high)];
        tb_dev_high = [tb_dev_high std(tb_vect_high)];

    end
end

%%
err_tb_mean_high = abs((tb_mean_high(end)-tb_mean_high(end-1))/tb_mean_high(end));
err_tb_dev_high = abs((tb_dev_high(end)-tb_dev_high(end-1))/tb_dev_high(end));

%%
figure
plot(tb_mean_low,'LineWidth',0.75)
hold on
plot(tb_mean_mid,'LineWidth',0.75)
plot(tb_mean_high,'LineWidth',0.75)
legend('Low pressure','Mid pressure', 'High pressure','Location','best')
title('Mean burning time')
xlabel('Iteration [-]')
ylabel('Time [s]')
grid on


figure
plot(tb_dev_low,'LineWidth',0.75)
hold on
plot(tb_dev_mid,'LineWidth',0.75)
plot(tb_dev_high,'LineWidth',0.75)
legend('Low pressure','Mid pressure', 'High pressure','Location','best')
title('Standard deviation of the burning time')
xlabel('Iteration [-]')
ylabel('Time [s]')
grid on

%%
rel_inc_tb_low = tb_dev_low(end)/tb_mean_low(end);
rel_inc_tb_mid = tb_dev_mid(end)/tb_mean_mid(end);
rel_inc_tb_high = tb_dev_high(end)/tb_mean_high(end);
rel_inc_a = Inc_a/a;
rel_inc_n = Inc_n/n;




