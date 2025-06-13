% 
% clear variables
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

min_expected_freq =0.3;
imageCount = 18000;
caseImages = 3600;
endpoint = '';
i = 5;
plot_type ='v';
modes = [1,2];
scaleFactor = 9.53;
    
        for run = 1:length(base)
            Calibrated = (fullfile(base{run}, 'CalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
            VelData = load(fullfile(Calibrated, sprintf('%05d.mat', 1)));
            [h, w] = size(VelData.piv_result(i).ux);
            
            mask = VelData.piv_result(i).b_mask;
            b_mask = mask;
            row_average = mean(mask, 2);
            index = find(row_average > 0.8, 1);  % separation row index
            % Preallocate cell arrays to store the structures
            numBatches = imageCount / caseImages;
            Data_above = cell(1, numBatches);
            Data_below = cell(1, numBatches);
            
            for batch = 1:numBatches
                dataloc_below = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                                         ['Cam' num2str(CameraNo)], 'Instantaneous', ...
                                         'Calibrated','Below', num2str(batch));
                dataloc_above = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                                         ['Cam' num2str(CameraNo)], 'Instantaneous', ...
                                         'Calibrated','Above', num2str(batch));
                Data_above{batch} = load(fullfile(dataloc_above, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
                Data_below{batch} = load(fullfile(dataloc_below, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
            end

            % Define locations for statistics and POD data
            statsloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                        ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
            dataloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                         ['Cam' num2str(CameraNo)], 'Instantaneous', ...
                         'Calibrated','Below');


            % (Optionally) adjust masks if needed
            mask_above = mask;
            mask_above(index:end, :) = 1;
            mask_below = mask;
            se = strel('disk', 3);  % Structuring element for dilation
            mask_below = imdilate(mask_below, se);  % Create buffer zone
            mask_below(1:index-1, :) = 1;
            
            % Load mean statistics (assumed to include mean_U and mean_V)
            StatsData = load(fullfile(statsloc, 'MeanStats16x16.mat'));
            mean_U = StatsData.mean_U;   % Expected to be vectorized for the region
            mean_V = StatsData.mean_V;
            u_prime_squared_calibrated = StatsData.U_prime_Uprime_mean;
            v_prime_squared_calibrated = StatsData.V_prime_Vprime_mean;
            TKE = 0.5 * (u_prime_squared_calibrated + v_prime_squared_calibrated);

            % Split TKE into the above and below regions
            TKE_masked_above = TKE(1:index, :);
            TKE_masked_below = TKE(index-1:end, :);
            Total_TKE_above = sum(TKE_masked_above(:));
            Total_TKE_below = sum(TKE_masked_below(:));

            % Compute uncertainty using the run-specific dt
            uncertainty = 0.1 / (scaleFactor/10^(-3) * dt_piv{run});
            b_mask_cavity = mask(index:end, :);
            Area_below = sum(~b_mask_cavity(:));
            Area_above = sum(~mask_above(:));

            Noise_TKE_below = Area_below * uncertainty^2;
            Noise_TKE_above = Area_above * uncertainty^2;

            % Calculate signal-to-noise ratios (SNR)
            SNR_below = Total_TKE_below / Noise_TKE_below;
            SNR_above = Total_TKE_above / Noise_TKE_above;

            % Determine energy percentages and modes required
            Noise_percentage_below = 100 * (1 / SNR_below);
            Noise_percentage_above = 100 * (1 / SNR_above);
            real_signal_below = 100 - Noise_percentage_below;
            real_signal_above = 100 - Noise_percentage_above;
            cumulativeEnergy_below = zeros(1, length(Data_below{1}.cumulativeEnergy));
            cumulativeEnergy_above = zeros(1, length(Data_below{1}.cumulativeEnergy));
            for batch = 1:(imageCount/caseImages)
                cumulativeEnergy_below = cumulativeEnergy_below + (Data_below{batch}.cumulativeEnergy);
                cumulativeEnergy_above = cumulativeEnergy_above + (Data_above{batch}.cumulativeEnergy);
            end
            cumulativeEnergy_below = cumulativeEnergy_below/5;
            cumulativeEnergy_above = cumulativeEnergy_above/5;
            k_below = find(cumulativeEnergy_below > real_signal_below, 1, 'first');
            k_above = find(cumulativeEnergy_above > real_signal_above, 1, 'first');
            
            disp(['The number of modes required to capture ', num2str(real_signal_below), '% of the energy is ', num2str(k_below), ' in cavity']);
            disp(['The number of modes required to capture ', num2str(real_signal_above), '% of the energy is ', num2str(k_above), ' above cavity']);
            
            
            
            dt = 1/sampling_rate(run);
            lower_top = -0.025; upper_top = 0.025;
            lower_bottom = -0.01; upper_bottom = 0.01;
            
            load(fullfile(dataloc, 'autocorr_results_split.mat'));
            
            
            
            xlim_l = -80;
            xlim_u = 80;
            ylim_u = 30;
            ylim_l = -45;
            for k = modes
                frame = ceil(results(k).period_lag/2);
                if isnan(results(k).period_lag)
    
                    continue
                end
                for p_t_loop = 1:2
                    if p_t_loop == 1
                        p_t = 'peak';
                    else
                        p_t = 'trough';
                    end

                    if strcmp(p_t, 'peak')
                        peak_data = load(fullfile(dataloc, ['peak_indices_Mode_split', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']));
                        indices = peak_data.peak_locations;
                         % Remove indices that would result in out-of-bound windows

                        % First, ensure the indices don't cross the overall data limits.
                        indices = indices((indices - frame >= 1) & (indices + frame <= imageCount));
                        
                        % Determine which set each index belongs to (assuming indices start at 1).
                        setIndices = ceil(indices / 3600);
                        
                        % Calculate the start and end indices for each set.
                        setStarts = (setIndices - 1) * 3600 + 1;
                        setEnds   = setIndices * 3600;
                        
                        % Only keep indices that, when subtracting/adding 'frame', remain within their set boundaries.
                        indices = indices((indices - frame >= setStarts) & (indices + frame <= setEnds));
                        
                        % Now, for each valid index, combine the indices from n-frame to n+frame.
                        combined_indices = [];
                        for n = indices
                            combined_indices = [combined_indices, (n - frame : n + frame)];
                        end

                    elseif strcmp(p_t, 'trough')
                        trough_data = load(fullfile(dataloc, ['trough_indices_Mode_split', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']));
                        
                        indices = trough_data.trough_locations;
                                                % First, ensure the indices don't cross the overall data limits.
                        indices = indices((indices - frame >= 1) & (indices + frame <= imageCount));
                        
                        % Determine which set each index belongs to (assuming indices start at 1).
                        setIndices = ceil(indices / 3600);
                        
                        % Calculate the start and end indices for each set.
                        setStarts = (setIndices - 1) * 3600 + 1;
                        setEnds   = setIndices * 3600;
                        
                        % Only keep indices that, when subtracting/adding 'frame', remain within their set boundaries.
                        indices = indices((indices - frame >= setStarts) & (indices + frame <= setEnds));
                        
                        % Now, for each valid index, combine the indices from n-frame to n+frame.
                        combined_indices = [];
                        for n = indices
                            combined_indices = [combined_indices, (n - frame : n + frame)];
                        end
                    end
                    if strcmp(plot_type, 'u')
                        mean_data =load(fullfile(base{run}, 'Statistics', num2str(imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', endpoint, ['MeanStats', num2str(16), 'x', num2str(16), '.mat']));
                        mean_field = mean_data.mean_U/ u_inf(run);
                    elseif strcmp(plot_type, 'v')
                        mean_data =load(fullfile(base{run}, 'Statistics', num2str(imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', endpoint,['MeanStats', num2str(16), 'x', num2str(16), '.mat']));
                        mean_field = mean_data.mean_V/ u_inf(run);
                    end
                    directory_path = fullfile(base{run}, 'CalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous');
                    Co_ords = load(fullfile(directory_path, "Co_ords"));
                    ycorners_top = [Co_ords.Co_ords(i).y(1, 1), Co_ords.Co_ords(i).y(index-1, end)];
                    xcorners_top = [Co_ords.Co_ords(i).x(1, 1), Co_ords.Co_ords(i).x(index-1, end)];
                    x = Co_ords.Co_ords(i).x;
                    y = Co_ords.Co_ords(i).y;


                    ycorners_bottom = [Co_ords.Co_ords(i).y(index, 1), Co_ords.Co_ords(i).y(end, end)];
                    xcorners_bottom = [Co_ords.Co_ords(i).x(index, 1), Co_ords.Co_ords(i).x(end, end)];
                    

                    numCells = 2* frame+1;
                    % This is preallocation for the shorter peak average video
                    developmentVideo_top = cell(1, numCells);
                    developmentVideo_bottom = cell(1, numCells);
                    developmentVideo_ux_top = cell(1, numCells);
                    developmentVideo_ux_bottom = cell(1, numCells);
                    developmentVideo_uy_top = cell(1, numCells);
                    developmentVideo_uy_bottom = cell(1, numCells);
                    developmentVideo_top_mean = cell(1, numCells);
                    developmentVideo_bottom_mean = cell(1, numCells);
                    developmentVideo_ux_top_mean = cell(1, numCells);
                    developmentVideo_uy_top_mean = cell(1, numCells);
                    developmentVideo_ux_bottom_mean = cell(1, numCells);
                    developmentVideo_uy_bottom_mean = cell(1, numCells);
                
                    % Initialize each cell with the relevant size array

                    for H = 1:numCells
                        developmentVideo_top{H} = zeros(size(VelData.piv_result(i).ux(1:index-1,:))); % Match Total_flow_top dimensions
                        developmentVideo_bottom{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_ux_top{H} = zeros(size(VelData.piv_result(i).ux(1:index-1,:))); % Match Total_flow_top dimensions
                        developmentVideo_uy_top{H} = zeros(size(VelData.piv_result(i).ux(1:index-1,:))); % Match Total_flow_top dimensions
                        developmentVideo_ux_bottom{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_uy_bottom{H} =zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_top_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_bottom_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_ux_top_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_uy_top_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_ux_bottom_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        developmentVideo_uy_bottom_mean{H} = zeros(size(VelData.piv_result(i).ux(index:end,:))); % Match Total_flow_bottom dimensions
                        
                    end
                  
                    count = 0;
                    energy_top = cell(numCells, length(indices));
                    energy_bottom = cell(numCells, length(indices));
                    counter =0;
                    index_count =1;
                    for idx = combined_indices
                        disp(idx)
                        batchNum = ceil(idx / 3600);
                        count = count +1;
                        imNo = idx - (batchNum - 1) * caseImages;
                        counter = counter +1;
                        if count > 2* frame+1
                            count =1;
                            index_count = index_count+1;
                        end
                        UV_rec_above = Data_above{batchNum}.U_svd(:, 1:k_above) * ...
                            (Data_above{batchNum}.S_svd(1:k_above, 1:k_above) * Data_above{batchNum}.V_svd(imNo, 1:k_above)');      


                        u_rec_above = UV_rec_above(1:h*w)+ Data_above{batchNum}.mean_u;
                        v_rec_above = UV_rec_above(h*w+1:end) + Data_above{batchNum}.mean_v;

                        u_rec_above = reshape(u_rec_above, [h, w]); % Size: (h x w)
                        v_rec_above = reshape(v_rec_above, [h, w]); % Size: (h x w)

                        u_rec_above_inst = UV_rec_above(1:h*w);
                        v_rec_above_inst = UV_rec_above(h*w+1:end);

                        u_rec_above_inst = reshape(u_rec_above_inst, [h, w]); % Size: (h x w)
                        v_rec_above_inst = reshape(v_rec_above_inst, [h, w]); % Size: (h x w)

                
                        %% Reconstruct the "below" (cavity) region
                        UV_rec_below = Data_below{batchNum}.U_svd(:, 1:k_below) * ...
                            (Data_below{batchNum}.S_svd(1:k_below, 1:k_below) * Data_below{batchNum}.V_svd(imNo, 1:k_below)');

                        u_rec_below = UV_rec_below(1:h*w)+ Data_below{batchNum}.mean_u;
                        v_rec_below = UV_rec_below(h*w+1:end) + Data_below{batchNum}.mean_v;

                        u_rec_below = reshape(u_rec_below, [h, w]); % Size: (h x w)
                        v_rec_below = reshape(v_rec_below, [h, w]); % Size: (h x w)

                        u_rec_below_inst = UV_rec_below(1:h*w);
                        v_rec_below_inst = UV_rec_below(h*w+1:end);

                        u_rec_below_inst = reshape(u_rec_below_inst, [h, w]); % Size: (h x w)
                        v_rec_below_inst = reshape(v_rec_below_inst, [h, w]); % Size: (h x w)


                
                        ux = [u_rec_above(1:index-1, :); u_rec_below(index:end, :)];
                        uy = [v_rec_above(1:index-1, :); v_rec_below(index:end, :)];

                        ux_inst = [u_rec_above_inst(1:index-1, :); u_rec_below_inst(index:end, :)];
                        uy_inst = [v_rec_above_inst(1:index-1, :); v_rec_below_inst(index:end, :)];

                        

                        ux = ux/ u_inf(run);
                        
                        uy = uy/ u_inf(run);
                        if strcmp(plot_type, 'u')
                            fluc = ux_inst/u_inf(run);
                        elseif strcmp(plot_type, 'v')
                            fluc = uy_inst/u_inf(run);
                        end
                        % total fluctating
                        fluc_top = fluc(1:index-1,:);
                        fluc_bottom = fluc(index:end,:);
                        ux_top = ux(1:index-1,:);
                        uy_top = uy(1:index-1,:);
                        ux_bottom = ux(index:end,:);
                        uy_bottom = uy(index:end,:);
                        energy_top{count,index_count} = sum(abs(fluc_top(:)));
                        energy_bottom{count,index_count} = sum(abs(fluc_bottom(:)));
                        developmentVideo_bottom{count} = developmentVideo_bottom{count} + fluc_bottom;
                        developmentVideo_top{count} = developmentVideo_top{count} + fluc_top;
                        developmentVideo_ux_top{count} = developmentVideo_ux_top{count} + ux_top;
                        developmentVideo_uy_top{count} = developmentVideo_uy_top{count} + uy_top;
                        developmentVideo_ux_bottom{count} = developmentVideo_ux_bottom{count} + ux_bottom;
                        developmentVideo_uy_bottom{count} = developmentVideo_uy_bottom{count} + uy_bottom;

                        % coeff = S_svd(k,k)* (V_svd(:,k));
                        % 
                        % time_coefficient = coeff(1:imageCount)';

                    end
                    count_no = 0;
                    for H = 1:numCells
                        count_no = count_no+1;
                        developmentVideo_top_mean{H} = developmentVideo_top{H} / length(indices);
                        developmentVideo_bottom_mean{H} = developmentVideo_bottom{H} / length(indices);
                        developmentVideo_ux_top_mean{H} = developmentVideo_ux_top{H} / length(indices);
                        developmentVideo_uy_top_mean{H} = developmentVideo_uy_top{H} / length(indices);
                        developmentVideo_ux_bottom_mean{H} = developmentVideo_ux_bottom{H} / length(indices);
                        developmentVideo_uy_bottom_mean{H} = developmentVideo_uy_bottom{H} / length(indices);
                    end
                    mp4_filename = fullfile(statsloc, ['FlowDevelopment-average', p_t, '-', plot_type, '-Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mp4']);
                    

                    v = VideoWriter(mp4_filename, 'MPEG-4');  % Create a VideoWriter object for MP4
                    v.FrameRate = 3;  % Set the frame rate (e.g., 10 frames per second)
                    open(v);  % Open the video writer
                    stat_path = (fullfile(statsloc,['energydevelopment', p_t, '-', plot_type, '-Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']));
                    save(stat_path, "energy_bottom","energy_top", "developmentVideo_top_mean","developmentVideo_bottom_mean","developmentVideo_ux_top_mean","developmentVideo_uy_top_mean","developmentVideo_uy_top_mean", "developmentVideo_ux_bottom_mean","developmentVideo_uy_bottom_mean")
                    for H = 1:numCells
                    
                        plotTimeInstanceGraph( b_mask, xcorners_top, ycorners_top, developmentVideo_top_mean{H}, lower_top, upper_top, ...
                            xcorners_bottom, ycorners_bottom, developmentVideo_bottom_mean{H}, lower_bottom, upper_bottom, ...
                            [], count_no, dt, xlim_l, xlim_u, ylim_l, ylim_u,'Fluctuating',index);
                        
                        % Capture frame
                        images = getframe(gcf);
                        writeVideo(v, images);

                        % Close figure
                        close(gcf);
                    end
                    close(v)
                 
                end

                

            end
        end
      