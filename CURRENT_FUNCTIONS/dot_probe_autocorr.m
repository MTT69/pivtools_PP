function dot_probe_autocorr(dot_probe_matfile, acquisition_frequency)
% DOT_PROBE_AUTOCORR - Plot symmetric autocorrelation of u and v probe time histories
% Inputs:
%   dot_probe_matfile - path to .mat file with u_time_history, v_time_history
%   acquisition_frequency - acquisition frequency in Hz

% Load data
data = load(dot_probe_matfile);
u_time_history = data.u_time_history;
v_time_history = data.v_time_history;
num_batches = data.num_batches;

% Find max lag (use length of shortest batch)
min_len = min(cellfun(@length, u_time_history));
max_lag = min_len - 1;
dt = 1/acquisition_frequency;
lags = (-max_lag:max_lag) * dt;

% Preallocate
u_acorr = zeros(num_batches, 2*max_lag+1);
v_acorr = zeros(num_batches, 2*max_lag+1);

% Compute symmetric autocorrelations for each batch
for b = 1:num_batches
    u = u_time_history{b};
    v = v_time_history{b};
    u = u(:) - mean(u); % zero mean
    v = v(:) - mean(v);
    u_corr = xcorr(u, max_lag, 'coeff');
    v_corr = xcorr(v, max_lag, 'coeff');
    u_acorr(b,:) = reshape(u_corr, 1, []);
    v_acorr(b,:) = reshape(v_corr, 1, []);
end

% Plot batch-wise symmetric autocorrelations
figure('Units','Normalized','OuterPosition',[0.1,0.1,0.6,0.8]);
for b = 1:num_batches
    subplot(num_batches,1,b);
    plot(lags, u_acorr(b,:), 'b-', 'LineWidth', 1.2); hold on;
    plot(lags, v_acorr(b,:), 'r-', 'LineWidth', 1.2);
    ylabel(sprintf('Batch %d', b));
    if b == 1
        title('Symmetric Autocorrelation of U (blue) and V (red) for Each Batch');
    end
    if b == num_batches
        xlabel('Time Lag (s)');
    end
    grid on;
    legend({'U','V'},'Location','northeast');
    xlim([min(lags), max(lags)]);
end

% Compute average symmetric autocorrelation
mean_u_acorr = mean(u_acorr,1);
mean_v_acorr = mean(v_acorr,1);

% Plot average symmetric autocorrelation
figure('Units','Normalized','OuterPosition',[0.2,0.2,0.5,0.5]);
plot(lags, mean_u_acorr, 'b-', 'LineWidth', 2); hold on;
plot(lags, mean_v_acorr, 'r-', 'LineWidth', 2);
xlabel('Time Lag (s)');
ylabel('Average Autocorrelation');
title('Average Symmetric Autocorrelation of U and V Across Batches');
legend({'U','V'},'Location','northeast');
grid on;
xlim([min(lags), max(lags)]);

end
