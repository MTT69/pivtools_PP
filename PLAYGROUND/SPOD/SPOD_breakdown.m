close all;
clear variables;

CameraNo = 1;
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
sampling_rate =[100,250,100,250,100,250,100,250,100,250,]; %hz 
u_inf = [0.34,1.18,0.34,1.18,0.34,1.18,0.34,1.18,0.34,1.18];
dt_piv = {3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6)}; % basic calibration dt - cell for different runs
imageCount = 18000;
caseImages = 3600;
i = 5;
plot_type ='v';
modes = [1,2];
scaleFactor = 9.53;
endpoint = "";
min_f = 0.3;





for run = 1:length(base)
    dataloc = fullfile(base{run}, 'CalibratedPIV', num2str(imageCount), ...
    ['Cam', num2str(CameraNo)], 'Instantaneous', endpoint);
    statsloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
        ['Cam' num2str(CameraNo)], 'Instantaneous', ...
        'Calibrated');

    VelData = load(fullfile(dataloc, sprintf('%05d.mat', 1)));
    Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
    ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
    xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
    [h, w] = size(VelData.piv_result(i).ux);
    mask = VelData.piv_result(i).b_mask;
    row_average = mean(mask, 2);
    index = find(row_average > 0.8, 1);
    se = strel('disk', 3); 
    
    % Apply dilation only to the rows starting from 'index' onwards
    mask(index:end, :) = imdilate(mask(index:end, :), se);  
    UV = zeros(imageCount,2*h,w);

    parfor ImNo = 1:imageCount  %% change
        % Load the data for each image
        VelData = load(fullfile(dataloc, sprintf('%05d.mat', ImNo)));
        u = VelData.piv_result(i).ux;
        v = VelData.piv_result(i).uy;
    
        % Apply the mask
        u(mask) = 0;
        v(mask) = 0;
    
        % Store the velocity fields
        UV(ImNo, :, :) = [u; v]; 
    
    end
    
    OPTS = struct();
    % Compute the mean field (temporal mean across frames)
    OPTS.mean = squeeze(mean(UV, 1));
    time_mean = OPTS.mean;
    segment_length = ceil(sampling_rate(run)/min_f);
    segment_length = min(caseImages, segment_length);
    nDFT= 2^nextpow2(segment_length); 
    overlap = nDFT/2;
    % Compute SPOD decomposition
    [L, P, f, ~, A] = spod(UV, nDFT, [], overlap, 1/sampling_rate(run),OPTS);
    % save spod stats in stats loc called ['SPOD_stats_', num2str(w), 'x', num2str(h), '.mat']
    save(fullfile(statsloc, ['SPOD_stats_', num2str(w), 'x', num2str(h), '.mat']), ...
    'L', 'P', 'f', 'A', 'nDFT', 'overlap', 'time_mean','-v7.3');



    
    % Figure for u velocity
    figure;
    set(gcf, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    count = 1;
    
    for fi = [3 4 5 6] % Iterate over frequency indices
        for mi = [1 ] % Iterate over mode indices
            mode_u = squeeze(P(fi, 1:h, :, mi)); % First half of height -> u
    
            % Subplot for u
            subplot(2, 2, count);
            contourf(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, real(mode_u), 11, 'EdgeColor', 'none');
            axis equal tight;
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 12);
            title(['$u: f=' num2str(f(fi),'%.2f') ',\ mode\ ' num2str(mi) '$'], ...
                'Interpreter', 'latex', 'FontSize', 12);
            colorbar;
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
            
            count = count + 1;
        end
    end
    
    colormap(jet);
    saveas(gcf, fullfile(statsloc, 'SPOD_decomposition_u.fig'));
    saveas(gcf, fullfile(statsloc, 'SPOD_decomposition_u.jpeg'));
    close gcf;
    
    % Figure for v velocity
    figure;
    set(gcf, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    count = 1;
    
    for fi = [1 2 3 4] % Iterate over frequency indices
        for mi = [1] % Iterate over mode indices
            mode_v = squeeze(P(fi, h+1:end, :, mi)); % Second half of height -> v
    
            % Subplot for v
            subplot(2, 2, count);
            contourf(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, real(mode_v), 11, 'EdgeColor', 'none');
            axis equal tight;
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
            ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 12);
            title(['$v: f=' num2str(f(fi),'%.2f') ',\ mode\ ' num2str(mi) '$'], ...
                'Interpreter', 'latex', 'FontSize', 12);
            colorbar;
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
            
            count = count + 1;
        end
    end
    
    colormap(jet);
    saveas(gcf, fullfile(statsloc, 'SPOD_decomposition_v.fig'));
    saveas(gcf, fullfile(statsloc, 'SPOD_decomposition_v.jpeg'));
    close gcf;




    figure
    loglog(f,L)
    xlabel('frequency'), ylabel('SPOD mode energy')
    
    % Compute cumulative modal energy for each frequency
    n = size(L, 1);  % Number of frequencies
    m = size(L, 2);  % Number of modes
    total_energy = sum(L(:));
    
    figure;
    hold on;
    for z = 1:10
        cumulative_energy = cumsum(L(z, :));
        Cumulative_energy_percentage = (cumulative_energy / total_energy) * 100;
        % Plot with a line and circle markers, increased line width and marker size
        plot(1:m, Cumulative_energy_percentage, '-o', ...
             'LineWidth', 1.5, 'MarkerSize', 6, ...
             'DisplayName', [num2str(f(z)), 'Hz']);
    end
    
    % Set axis labels and title with LaTeX interpreter and custom font size
    xlabel('Mode', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Cumulative Modal Energy Percentage', 'Interpreter', 'latex', 'FontSize', 12);
    title('Cumulative Modal Energy for Each Frequency', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Display the legend with LaTeX formatting
    legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'best');
    
    % Enhance the plot appearance
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    hold off;

    %save as .fig and .jpeg
    saveas(gcf, fullfile(statsloc, 'Cumulative_modal_energy.fig'));
    saveas(gcf, fullfile(statsloc, 'Cumulative_modal_energy.jpeg'));
    close gcf;
end



   



