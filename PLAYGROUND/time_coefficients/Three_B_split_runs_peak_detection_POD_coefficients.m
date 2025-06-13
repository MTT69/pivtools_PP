
    
base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};
samplingRate =[100,250,100,250,100,250,100,250,100,250]; %hz 
min_expected_freq =0.3;
imageCount = 18000;
caseImages = 3600;
endpoint = '';
x_loc = -5;
scaleFactor = 9.53;
CameraNo =1;
close all
i =5;


for run = 1:length(base)


    
    uncalibrated = (fullfile(base{run}, 'UncalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
    VelData = load(fullfile(uncalibrated, (sprintf('%05d.mat', 1))));
    [h,w] = size(VelData.piv_result(i).ux);
    timecoefficient = zeros(3,imageCount);
    for batch = 1 : (imageCount/caseImages)
        dataloc = (fullfile(base{run}, 'Statistics', num2str(imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated','Below',num2str(batch)));
        Data = load(fullfile(dataloc, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
        
        for k = 1:3
            mode = k;
            time_coefficient = Data.S_svd(mode,mode)* (Data.V_svd(:,mode)');
            time_coefficient = time_coefficient - mean(time_coefficient);
            timecoefficient(k,((batch-1)*caseImages+1):batch*caseImages) = time_coefficient;
        end
    end


    for k = 1:3
        mode = k;
        dt = 1/samplingRate(run);
        num_frames = imageCount;
        A_prime_detrended = timecoefficient(k,:);
        time = (0:imageCount-1) * dt; % Time array (in seconds)
        [peaks, peak_locations] = findpeaks(A_prime_detrended, 'MinPeakWidth', 40);

        % Find troughs (invert the signal to treat troughs as peaks)
        [troughs, trough_locations] = findpeaks(-A_prime_detrended, 'MinPeakWidth', 40);

        % Plot the time coefficient against time
        figure;
        plot(time, A_prime_detrended, 'b', 'LineWidth', 1.5);
        hold on;

        % Add vertical red lines for acquisition points
        for z = 1:(imageCount/caseImages - 1)
            xline(z * caseImages * dt, 'r--', 'LineWidth', 1.2);
        end

        % Mark peaks in blue
        plot(time(peak_locations), A_prime_detrended(peak_locations), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Peaks');

        % Mark troughs in red
        plot(time(trough_locations), A_prime_detrended(trough_locations), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Troughs');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        % Customize the plot
        xlabel('Time (s)');
        ylabel('Time Coefficient');
        title('Time Coefficient vs Time with Peaks and Troughs');
        grid on;

        % Add legend
        hold off;
        dataloc = (fullfile(base{run}, 'Statistics', num2str(imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated','Below'));

        % save the figure file name needs to contain mode number as well as h and w
        filename = ['TimeCoefficient_Mode_split', num2str(mode), '_', num2str(w), 'x', num2str(h)]; 
        fullpath_jpg = fullfile(dataloc, [filename, '.jpg']);
        fullpath_epsc = fullfile(dataloc, [filename, '.epsc']);
        fullpath_fig = fullfile(dataloc, [filename, '.fig']);
        saveas(gcf, fullpath_jpg, 'jpg');
        saveas(gcf, fullpath_epsc, 'epsc');
        saveas(gcf, fullpath_fig, 'fig');
        close(gcf);
        % also save the peak and trough indices again with mode number and h and w
        save(fullfile(dataloc, ['trough_indices_Mode_split', num2str(mode), '_', num2str(w), 'x', num2str(h), '.mat']), 'troughs', 'trough_locations')
        save(fullfile(dataloc, ['peak_indices_Mode_split', num2str(mode), '_', num2str(w), 'x', num2str(h), '.mat']),  'peaks', 'peak_locations')

 
    end

end







