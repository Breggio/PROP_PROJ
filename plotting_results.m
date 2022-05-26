clear
load("MC_results.mat");
figure(1);

% Plotting cumulative mean
subplot(1,2,1);
plot(Cumulative_mean_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex', 'FontSize',21);
ylabel('Cumulative mean [N]', 'Interpreter', 'latex', 'FontSize',21);
title('\textbf{Cumulative mean of thurst}', 'Interpreter', 'latex', 'FontSize',21);

% Plotting cumulative standard deviation
subplot(1,2,2);
plot(Cumulative_standard_deviation_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex', 'FontSize',21);
ylabel('Cumulative standard deviation [N]', 'Interpreter', 'latex', 'FontSize',21);
title('\textbf{Cumulative standard deviation of thrust}', 'Interpreter', 'latex', 'FontSize',21);

sgtitle('\textbf{Monte Carlo analysis}', 'Interpreter', 'latex');
