% Assuming L is your energy array with size n x m (n frequencies, m modes)
% f is the frequency array corresponding to the rows of L.
close all
n = size(L, 1);  % Number of frequencies
m = size(L, 2);  % Number of modes

% --- Plot 1: Cumulative total energy percentage of the frequencies ---

% Compute total energy for each frequency (sum over modes)
frequency_energy = sum(L, 2);

total_energy = sum(frequency_energy);

relative_energy = frequency_energy /total_energy *100;

% Plot cumulative energy percentage vs frequencies
figure;
plot(f, frequency_energy, '-o', 'LineWidth', 2);
hold on;



grid on;
xlabel('Frequency hz');
ylabel('Total Energy fraction (%)');
title('Cumulative Energy Percentage of Frequencies');

% % Save the plot
% saveas(gcf, fullfile(directory, 'Cumulative_Frequency_Energy.jpg'));
% saveas(gcf, fullfile(directory, 'Cumulative_Frequency_Energy.epsc'));
% saveas(gcf, fullfile(directory, 'Cumulative_Frequency_Energy.fig'));
% close all;

% --- Plot 2: Cumulative modal energy for the top 10 most energetic frequencies ---

% Sort frequencies by their total energy and select top 10
[~, sorted_indices] = sort(frequency_energy, 'descend');
top_frequencies = sorted_indices(1:min(10, n));

% Initialize a figure for the modal energy plots
figure;
hold on;

for idx = 1:length(top_frequencies)
    freq_idx = top_frequencies(idx);

    % Compute cumulative modal energy for this frequency
    cumulative_modal_energy = cumsum(L(freq_idx, :));

    % Plot the cumulative modal energy
    plot(1:m, cumulative_modal_energy, 'DisplayName', sprintf('f = %.2f', f(freq_idx)));
end

grid on;
xlabel('Mode');
ylabel('Cumulative Modal Energy');
title('Cumulative Modal Energy for Top 10 Frequencies');
legend show;

% % Save the plot
% saveas(gcf, fullfile(directory, 'Top10_Cumulative_Modal_Energy.jpg'));
% saveas(gcf, fullfile(directory, 'Top10_Cumulative_Modal_Energy.epsc'));
% saveas(gcf, fullfile(directory, 'Top10_Cumulative_Modal_Energy.fig'));
% close all;