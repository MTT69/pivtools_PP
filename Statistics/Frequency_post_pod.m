function Frequency_post_pod(setup,Type,CameraNo,endpoint,domain, sampling_rate, frame_range,modes,min_expected_freq)
    if setup.pipeline.POST_POD_stats
        
        for i = setup.instantaneous.runs

        
            dataloc = (fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', Type, endpoint, domain));
            uncalibrated = (fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
            VelData = load(fullfile(uncalibrated, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
            [h,w] = size(VelData.piv_result(i).ux);
            Data = load(fullfile(dataloc, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
            A = Data. A;
            clear Data

            for k = modes
                mode = modes(k);
                dt = 1/sampling_rate;
                num_frames = setup.imProperties.imageCount;

                time_coefficient = A(1:setup.imProperties.imageCount, mode);
                time = (0:setup.imProperties.imageCount-1) * dt; % Time array (in seconds)
                % Detrend the data for peak/trough detection
                A_prime_detrended = time_coefficient - mean(time_coefficient);
                [peaks, peak_locations] = findpeaks(A_prime_detrended, 'MinPeakWidth', 40);

                % Find troughs (invert the signal to treat troughs as peaks)
                [troughs, trough_locations] = findpeaks(-A_prime_detrended, 'MinPeakWidth', 40);

                % Extended range around peaks
                extended_peak_locations = extendLocations(peak_locations, frame_range, num_frames);

                % Extended range around troughs
                extended_trough_locations = extendLocations(trough_locations, frame_range, num_frames);

                % Plot the time coefficient against time
                figure;
                plot(time, time_coefficient, 'b', 'LineWidth', 1.5);
                hold on;

                % Add vertical red lines for acquisition points
                for z = 1:(setup.imProperties.imageCount/setup.imProperties.caseImages - 1)
                    xline(z * setup.imProperties.caseImages * dt, 'r--', 'LineWidth', 1.2);
                end

                % Mark peaks in blue
                plot(time(peak_locations), time_coefficient(peak_locations), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Peaks');

                % Mark troughs in red
                plot(time(trough_locations), time_coefficient(trough_locations), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Troughs');
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                % Customize the plot
                xlabel('Time (s)');
                ylabel('Time Coefficient');
                title('Time Coefficient vs Time with Peaks and Troughs');
                grid on;

                % Add legend
                hold off;
                % save the figure file name needs to contain mode number as well as h and w
                filename = ['TimeCoefficient_Mode', num2str(mode), '_', num2str(w), 'x', num2str(h)]; 
                fullpath_jpg = fullfile(dataloc, [filename, '.jpg']);
                fullpath_epsc = fullfile(dataloc, [filename, '.epsc']);
                fullpath_fig = fullfile(dataloc, [filename, '.fig']);
                saveas(gcf, fullpath_jpg, 'jpg');
                saveas(gcf, fullpath_epsc, 'epsc');
                saveas(gcf, fullpath_fig, 'fig');
                close(gcf);
                % also save the peak and trough indices again with mode number and h and w
                save(fullfile(dataloc, ['trough_indices_Mode', num2str(mode), '_', num2str(w), 'x', num2str(h), '.mat']), 'troughs', 'trough_locations', "extended_trough_locations")
                save(fullfile(dataloc, ['peak_indices_Mode', num2str(mode), '_', num2str(w), 'x', num2str(h), '.mat']),  'peaks', 'peak_locations', "extended_peak_locations")

                %% Welch's Method for PSD Estimation
                segment_length = ceil(sampling_rate/min_expected_freq);
                segment_length = min(setup.imProperties.caseImages, segment_length);
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
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                grid on;

                % Find the predominant frequency
                [max_power, idx] = max(pxx); % Find the peak in the PSD
                predominant_frequency = f(idx);
                disp(['Mode: ', num2str(mode), ', Run: ', num2str(i), ', Predominant Frequency: ', num2str(predominant_frequency), ' Hz']);
                % save the figure file name needs to contain mode number as well as h and w
                filename = ['PSD_Mode', num2str(mode), '_', num2str(w), 'x', num2str(h)];
                fullpath_jpg = fullfile(dataloc, [filename, '.jpg']);
                fullpath_epsc = fullfile(dataloc, [filename, '.epsc']);
                fullpath_fig = fullfile(dataloc, [filename, '.fig']);
                saveas(gcf, fullpath_jpg, 'jpg');
                saveas(gcf, fullpath_epsc, 'epsc');
                saveas(gcf, fullpath_fig, 'fig');
                close(gcf);
                % also save the PSD data
                save(fullfile(dataloc, ['PSD_Mode', num2str(mode), '_', num2str(w), 'x', num2str(h), '.mat']), 'pxx', 'f', 'predominant_frequency')
            end
        end
    end
    function [extended_locations] = extendLocations(locations, frame_range, num_frames)
        extended_locations = cell(length(locations), 1);
        for j = 1:length(locations)
            index = locations(j);
            start_index = max(1, index - frame_range);
            end_index = min(num_frames, index + frame_range);
            extended_locations{j} = start_index:end_index;
        end
        extended_locations = unique([extended_locations{:}]);  % Flatten and remove duplicates
    end
end






