load("results_of_monte_carlo.mat");
figure(1);

% Plotting cumulative mean
subplot(1,2,1);
plot(Cumulative_mean_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative mean [N]', 'Interpreter', 'latex');
title('\textbf{Cumulative mean of thurst}', 'Interpreter', 'latex');

% Plotting cumulative standard deviation
subplot(1,2,2);
plot(Cumulative_standard_deviation_thrust, 'LineWidth', 1.5);
hold on
grid on
xlabel('Monte Carlo iterations', 'Interpreter', 'latex');
ylabel('Cumulative standard deviation [N]', 'Interpreter', 'latex');
title('\textbf{Cumulative standard deviation of thrust}', 'Interpreter', 'latex');
