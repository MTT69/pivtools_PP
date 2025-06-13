% Check if data 'A' exists, if not, load the specified file
try
    A;
catch
    load('E:\Processed_PIV\90degree_400light_100hz_3000dt\Statistics\18000\Cam1\Instantaneous\Calibrated\RPCA\Below\POD_stats_319x149.mat');
end

% Clear unnecessary variables
clear A_s LAM_s PHI PHI_s lambda_s ilam_s
mode = 1;
% Define parameters
sampling_rate = 100; % Hz
t = 18000;           % Total points
run_count = 3600;    % Points per run
dt = 1 / sampling_rate; % Time step (seconds per sample)
frame_range = 20;    % Frame range around peaks/troughs
num_frames = t;      % Total frames in the data

% Extract required data
time_coefficient = A(1:t, mode);

% Create a time vector
time = (0:t-1) * dt; % Time array (in seconds)

% Detrend the data for peak/trough detection
A_prime_detrended = time_coefficient - mean(time_coefficient);

% Find peaks
[peaks, peak_locations] = findpeaks(A_prime_detrended, 'MinPeakWidth', 40);

% Find troughs (invert the signal to treat troughs as peaks)
[troughs, trough_locations] = findpeaks(-A_prime_detrended, 'MinPeakWidth', 40);

% Extended range around peaks
extended_peak_locations = cell(length(peak_locations), 1);
for j = 1:length(peak_locations)
    peak_index = peak_locations(j);
    start_index = max(1, peak_index - frame_range);
    end_index = min(num_frames, peak_index + frame_range);
    extended_peak_locations{j} = start_index:end_index;
end
extended_peak_locations = unique([extended_peak_locations{:}]);  % Flatten and remove duplicates
% Extended range around troughs
extended_trough_locations = cell(length(trough_locations), 1);
for j = 1:length(trough_locations)
    trough_index = trough_locations(j);
    start_index = max(1, trough_index - frame_range);
    end_index = min(num_frames, trough_index + frame_range);
    extended_trough_locations{j} = start_index:end_index;
end
extended_trough_locations = unique([extended_trough_locations{:}]);  % Flatten and remove duplicates

% Plot the time coefficient against time
figure;
plot(time, time_coefficient, 'b', 'LineWidth', 1.5);
hold on;

% Add vertical red lines for acquisition points
for i = 1:(t/run_count - 1)
    xline(i * run_count * dt, 'r--', 'LineWidth', 1.2);
end

% Mark peaks in blue
plot(time(peak_locations), time_coefficient(peak_locations), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Peaks');

% Mark troughs in red
plot(time(trough_locations), time_coefficient(trough_locations), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Troughs');

% Customize the plot
xlabel('Time (s)');
ylabel('Time Coefficient');
title('Time Coefficient vs Time with Peaks and Troughs');
grid on;

% Add legend
legend('Time Coefficient', 'Run Start', 'Peaks', 'Troughs');
hold off;

%% Welch's Method for PSD Estimation
% Define parameters for Welch's method
segment_length = run_count;  % 3600 samples per segment
overlap = segment_length / 2; % 50% overlap
nfft = 2^nextpow2(segment_length); % Number of FFT points

% Compute the PSD using Welch's method
[pxx, f] = pwelch(time_coefficient, hamming(segment_length), overlap, nfft, sampling_rate);

% Plot the PSD
figure;
plot(f, 10*log10(pxx), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density using Welch Method');
grid on;

% Find the predominant frequency
[max_power, idx] = max(pxx); % Find the peak in the PSD
predominant_frequency = f(idx);
disp(['Predominant Frequency: ', num2str(predominant_frequency), ' Hz']);


save(fullfile(pwd,'trough_indices.mat'), 'extended_trough_locations')
save(fullfile(pwd,'peak_indices.mat'),  'extended_peak_locations')



