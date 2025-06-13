function frequencies = process_and_plot(data, dataLabel, save_dir, samplingRate, segment_length, overlap, nfft)
    % Function to compute and plot time series, autocorrelation, and PSD
    % data: The time series data to be processed
    % dataLabel: The name used in plot titles and filenames
    % save_dir: Directory to save the plots
    % samplingRate, segment_length, overlap, nfft: Parameters for PSD computation
    frequencies =nan;
    % Create time series plot
    nanRatio = sum(isnan(data)) / numel(data);
    if nanRatio > 0.3
        disp(['Too many NaNs for ' dataLabel,' NaN fraction: ' [num2str(nanRatio)]]);
        return
    end
    
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plot(1:length(data), data, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    title([dataLabel, ' vs Time']);
    xlabel('Time Step');
    ylabel(dataLabel);
    xlim([0, 500]);

    % Save figure
    saveas(gcf, fullfile(save_dir, [dataLabel, '.fig']));
    saveas(gcf, fullfile(save_dir, [dataLabel, '.jpg']));
    close(gcf);

    xc_length = 2 * 3600 - 1;
    
    % Initialize an accumulator for the autocorrelation and a counter.
    local_ac_sum = zeros(1, xc_length);
    
    for seg = 1:5
        % Define index range for this segment.
        idx_start = (seg - 1) * 3600 + 1;
        idx_end   = seg * 3600;
        
        % Extract the segment from the data.
        seg_data = data(idx_start:idx_end)';
        
        % Preprocess the segment:
        % 1. Fill missing data.
        % 2. Remove the mean.
        seg_data = fillmissing(seg_data, 'linear');
        seg_data = seg_data - mean(seg_data);
        
        % Compute autocorrelation for this segment.
        ac_seg = xcorr(seg_data, 'unbiased');
        
        % Accumulate the autocorrelation.
        local_ac_sum = local_ac_sum + ac_seg;
    end
    
    % Average the autocorrelation over all segments.
    autocorr_data = local_ac_sum / 5;

    % Plot autocorrelation
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plot(autocorr_data);
    title(['Autocorrelation of ', dataLabel]);
    xlabel('Lag');
    ylabel('Autocorrelation');
    grid on;
    xlim([3100,4100]);

    % Save autocorrelation plot
    saveas(gcf, fullfile(save_dir, ['Autocorrelation_', dataLabel, '.fig']));
    saveas(gcf, fullfile(save_dir, ['Autocorrelation_', dataLabel, '.jpg']));
    close(gcf);

    % Compute Power Spectral Density (PSD)
    data = fillmissing(data, 'linear'); % Handle missing data if any
    data = data - mean(data); % Remove the mean
    
    [pxx, f] = pwelch(data, hamming(segment_length), overlap, nfft, samplingRate);

    % Plot PSD
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plot(f, 10 * log10(pxx), 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(['Power Spectral Density (PSD) of ', dataLabel]);
    grid on;
    xlim([0, 5]);

    % Save PSD plot
    saveas(gcf, fullfile(save_dir, ['PSD_', dataLabel, '.fig']));
    saveas(gcf, fullfile(save_dir, ['PSD_', dataLabel, '.jpg']));
    close(gcf);

    % Find the frequency corresponding to the maximum energy
    [~, idx] = max(pxx);
    frequencies = f(idx); % Frequency at maximum power
end