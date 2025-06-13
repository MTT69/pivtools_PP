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
samplingRate =[100,250,100,250,100,250,100,250,100,250,]; %hz 
dt = {3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),}; % basic calibration dt - cell for different runs
min_expected_freq =0.3;
imageCount = 18000;
caseImages = 3600;
endpoint = '';
run = 5;
x_loc = [-5,-100];
scaleFactor = 9.53;
CameraNo =1;
close all


for i = 1:length(base)
    disp(['Processing base directory: ', base{i}]);
    disp(datetime('now', 'Format', 'yyyy, MM, dd'))
    statsloc = fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    % Load Co_ords.mat
    VelData = load(fullfile(statsloc, 'MeanStats16x16.mat'));
    u_prime_squared_calibrated = VelData.U_prime_Uprime_mean;
    v_prime_squared_calibrated = VelData.V_prime_Vprime_mean;
    mean_U = VelData.mean_U;
    mean_V = VelData.mean_V;


    TKE = 0.5*(u_prime_squared_calibrated + v_prime_squared_calibrated);
    POD_stats = load(fullfile(statsloc, 'Below', 'POD_stats_319x149.mat')); % harded coded oops
    
    h=POD_stats.h;
    w=POD_stats.w;
    dataloc = (fullfile(base{i}, 'CalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
    % Load Co_ords.mat
    VelData = load(fullfile(dataloc, [num2str(sprintf('%05d.mat', 1))]));
    Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
    b_mask = VelData.piv_result(run).b_mask;
    se = strel('disk', 3);  % Structuring element for dilation
    b_mask_dilated = imdilate(b_mask, se);  % Create buffer zone
    row_means = mean(b_mask,2);
    first_row_index = find(row_means > 0.8, 1, 'first');
    y = Co_ords.Co_ords(run).y;
    x = Co_ords.Co_ords(run).x;
    ycorners = [Co_ords.Co_ords(run).y(1,1), Co_ords.Co_ords(run).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
    xcorners = [Co_ords.Co_ords(run).x(1,1), Co_ords.Co_ords(run).x(end, end)];

    TKE_masked = TKE;
    TKE_masked = TKE_masked(first_row_index:end,:);
    Total_TKE = sum(TKE_masked(:));
    % we are working in a calibrated space to the 0.1 error needs to be in relation to m/s
    uncertainty = 0.1/(scaleFactor/10^(-3)*dt{i});
    b_mask_cavity = b_mask(first_row_index:end,:);
    cavity_mask = (~b_mask_cavity);
    Area = sum(~b_mask_cavity(:));
    Noise_TKE = (Area*uncertainty^2);
    % Calculate the SNR
    SNR = Total_TKE/ Noise_TKE;
    % Noise percentage
    Noise_percentage = 100*(1/SNR);

    real_signal = 100- Noise_percentage;

    % Find the index of the first value that is greater than the real signal
    k = find(POD_stats.cumulativeEnergy > real_signal, 1, 'first');
    disp(['The number of modes required to capture ', num2str(real_signal), '% of the energy is ', num2str(k)]);
    POD_rebuild = struct();
    POD_rebuild.modes_required = k;
    POD_rebuild.SNR = SNR;
    POD_rebuild.signal = real_signal;


    vortexCentersX_all = nan(imageCount, 2);
    vortexCentersY_all = nan(imageCount, 2);
    rotation_all = nan(imageCount, 2);
    v_flux_all = nan(imageCount, 8);
    deltaStar_all = nan(imageCount, 2);
    vortexSize_all = nan(imageCount, 2);
    momentumThickness_all = nan(imageCount, 2);
    convective_velocity = nan(imageCount, 8);


    for imNo = 1: imageCount
        %% recreate flow 
        dataloc = (fullfile(base{i}, 'CalibratedPIV', num2str(imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
        % Load Co_ords.mat
        VelData = load(fullfile(dataloc, [num2str(sprintf('%05d.mat', imNo))]));

        u = VelData.piv_result(run).ux;
        v = VelData.piv_result(run).uy;

        UV_rec = POD_stats.U_svd(:, 1:k) * ((POD_stats.S_svd(1:k, 1:k)) * POD_stats.V_svd(imNo, 1:k)'); 
%         S_svd(mode,mode)* (V_svd(:,mode));

        % Split the reconstructed data into u and v components ( this is
        % instantaneous fluctuation)
        u_rec = UV_rec(1:h*w) + POD_stats.mean_u;       % First half corresponds to u
        v_rec = UV_rec(h*w+1:end) + POD_stats.mean_v;   % Second half corresponds to v

        u_rec = reshape(u_rec, [h, w]); % Size: (h x w)
        v_rec = reshape(v_rec, [h, w]); % Size: (h x w)
        
        u(first_row_index+1:end,:) = u_rec(first_row_index+1:end,:);
        v(first_row_index+1:end,:) = v_rec(first_row_index+1:end,:);



        %% Boundary layer 
        % Select only y > 0 (upper boundary layer)
        positive_y_indices = y(:,1) > 0;

        % Define the freestream velocity (99th percentile of the x-velocity)
        percentile_speed = prctile(u(positive_y_indices, :), 99, 'all');

        deltaStar_temp = zeros(1, length(x_loc));

        momentumThickness_temp = zeros(1, length(x_loc));

        for j = 1:length(x_loc)
            % Find the closest x-column to the desired x-location
            [~, x_index] = min(abs(x(1,:) - x_loc(j)));
            
            % Define a range of columns around this x-location to smooth the profile
            col_range = max(1, x_index-3) : min(size(x, 2), x_index+3);
            
            % Extract and normalize the velocity profile (for y > 0)
            velocity_profile = mean(u(positive_y_indices, col_range), 2) / percentile_speed;
            
            % Compute the displacement thickness δ* using trapezoidal integration
            deltaStar_all(imNo,j) = trapz(flipud(y(positive_y_indices, 1)), (1 - velocity_profile));
            
            % Compute the momentum thickness θ using trapezoidal integration
            momentumThickness_all(imNo,j) = trapz(flipud(y(positive_y_indices, 1)), velocity_profile .* (1 - velocity_profile));
        end

        %% gammaField
        gammaField = gamma1(x, y, u, v, 10);
        gammaField(b_mask) = 0;
        gammaField(1:20,:)=0;
        gammaFieldLocation = gammaField;
        % Set the values where b_mask is true to NaN


        gammaFieldLocation(b_mask) = NaN;
        % Threshold for vortex detection
        threshold = 0.5;
        
        % Identify Positive (Counterclockwise) Vortex Centers
        BW_pos = imregionalmax(gammaField); % Local maxima
        BW_pos = BW_pos & ~b_mask_dilated & (y <= 0);
        
        
        
        % Identify Negative (Clockwise) Vortex Centers
        BW_neg = imregionalmin(gammaField); % Local minima
        BW_neg = BW_neg & ~b_mask_dilated & (y <= 0);
        
        % Extract vortex center coordinates and gamma values
        vortexCentersX_pos = x(BW_pos)';
        vortexCentersY_pos = y(BW_pos)';
        gamma2Values_pos = gammaField(BW_pos)';
        
        vortexCentersX_neg = x(BW_neg)';
        vortexCentersY_neg = y(BW_neg)';
        gamma2Values_neg = gammaField(BW_neg)';

        if isempty(vortexCentersX_pos) % No positive vortex found
            vortexPosX = NaN;
            vortexPosY = NaN;
        else
            [~, idx_max_pos] = max(gamma2Values_pos); % Strongest positive vortex
            vortexPosX = vortexCentersX_pos(idx_max_pos);
            vortexPosY = vortexCentersY_pos(idx_max_pos);
        end
        
        if isempty(vortexCentersX_neg) % No negative vortex found
            vortexNegX = NaN;
            vortexNegY = NaN;
        else
            [~, idx_max_neg] = min(gamma2Values_neg); % Strongest negative vortex
            vortexNegX = vortexCentersX_neg(idx_max_neg);
            vortexNegY = vortexCentersY_neg(idx_max_neg);
        end
        
        % Vortex centres ' maximum value' 
        vortexCentersX = [vortexPosX, vortexNegX];
        vortexCentersY = [vortexPosY, vortexNegY];
        for k = 1:2
            % Skip if invalid or missing
           if isnan(vortexCentersX(k))
                % Assign NaN only to the missing vortex, but continue processing the other
                vortexCentersX_all(imNo,k) = NaN;
                vortexCentersY_all(imNo,k) = NaN;
                rotation_all(imNo,k) = NaN;
                vortexSize_all(imNo,k) = NaN;
                continue; % Move to the next vortex without breaking the loop
           end
      

            % 1) Find the closest grid point for this vortex center
            [~, centerIdx] = min((x(:) - vortexCentersX(k)).^2 + (y(:) - vortexCentersY(k)).^2);
            peakGamma = gammaField(centerIdx); 
            contourLevel = abs(0.5* peakGamma); % 2/pi 
    
            % 2) Create a binary mask for this vortex region
            regionMask = (abs(gammaField) >= contourLevel) & (abs(gammaField) >= 0.50);
    
            regionMask = bwareaopen(regionMask,30); % Remove small fragmented components
            % Create a kernel (5x1 column filter) to check the vertical neighborhood
            kernel = [1; 1; 1; 1; 1];  % 5x1 kernel, sliding over the vertical direction
            
            % Apply convolution to count the number of non-zero neighbors
            countNonZeroNeighbors = conv2(regionMask == 0, kernel, 'same'); 
            
            % Create a mask where 3 or more neighbors are zero (i.e., the count of zero neighbors is >= 3)
            invalidMask = countNonZeroNeighbors >= 3;
            
            % Apply the invalid mask to set those rows in regionalMask to zero
            regionMask(invalidMask) = 0;

            
            CC = bwconncomp(regionMask); % Get connected components
            numPixels = cellfun(@numel, CC.PixelIdxList); % Count number of pixels in each component
            
            % Find the indices of the two largest components
            [~, sortedIndices] = maxk(numPixels, 2); 
            
            % Create a new binary mask with only the two largest components
            newRegionMask = false(size(regionMask));
            for L = 1:length(sortedIndices)
                newRegionMask(CC.PixelIdxList{sortedIndices(L)}) = true;
            end
            
            % Multiply by gammaField to retain original values
            regionMask = newRegionMask .* gammaField;
    
    
            
            % 3) Identify which connected component contains the vortex
            % center for this index
            CC = bwconncomp(regionMask);
            componentIndex = [];
            for compNum = 1:CC.NumObjects
                if ismember(centerIdx, CC.PixelIdxList{compNum})
                    componentIndex = compNum;
                    break;
                end
            end
            if isempty(componentIndex)
                vortexCentersX_all(imNo,k) = NaN;
                vortexCentersY_all(imNo,k) = NaN;
                rotation_all(imNo,k) = NaN;
                vortexSize_all(imNo,k) = NaN;
                continue; % Move to the next vortex without breaking the loop
            end
            isolatedMask = zeros(size(regionMask));  % Initialize empty mask
            isolatedMask_weights = isolatedMask;
            vorticityMask = zeros(size(regionMask));
            pixelIndices = CC.PixelIdxList{componentIndex};  % Get the pixel indices for the identified component
            
            
            isolatedMask_weights(pixelIndices) = gammaField(pixelIndices);  % Keep only the identified component
            isolatedMask = abs(isolatedMask_weights) >0;
            [dx,~] = gradient(x);
            [~,dy] = gradient(y);
            dy = -dy;  % because y was flipped
            pixel_area = (dx(1,1) * dy(1,1))*1*(10^-6);
            
       
              
            % 5) Calculate the centroid
            regionStats = regionprops(isolatedMask, isolatedMask_weights, 'WeightedCentroid');
            ctd = regionStats.WeightedCentroid;
            % Extract the fractional part and integer part of the centroid
            rowLower = floor(ctd(2));   % The integer part of the row (y-coordinate)
            colLower = floor(ctd(1));   % The integer part of the column (x-coordinate)
            
            rowUpper = rowLower + 1;    % The next row (y-coordinate)
            colUpper = colLower + 1;    % The next column (x-coordinate)
            
            fractionalRow = ctd(2) - rowLower;  % The fractional part of the row
            fractionalCol = ctd(1) - colLower;  % The fractional part of the column
            
            % Interpolate the real-world coordinates (x and y) at the centroid location
            % Assuming x and y are the real-world calibrated meshgrids
            
            % Interpolate x-coordinate
            vortexCentersX_all(imNo,k) = x(rowLower, colLower) + fractionalRow * (x(rowUpper, colLower) - x(rowLower, colLower)) + ...
                            fractionalCol * (x(rowLower, colUpper) - x(rowLower, colLower)) + ...
                            fractionalRow * fractionalCol * (x(rowUpper, colUpper) - x(rowUpper, colLower) - x(rowLower, colUpper) + x(rowLower, colLower));
            


            
            % Interpolate y-coordinate
            vortexCentersY_all(imNo,k) = y(rowLower, colLower) + fractionalRow * (y(rowUpper, colLower) - y(rowLower, colLower)) + ...
                            fractionalCol * (y(rowLower, colUpper) - y(rowLower, colLower)) + ...
                            fractionalRow * fractionalCol * (y(rowUpper, colUpper) - y(rowUpper, colLower) - y(rowLower, colUpper) + y(rowLower, colLower));
            if sqrt((vortexCentersY_all(imNo,k)-vortexCentersY(k))^2 + (vortexCentersX_all(imNo,k)-vortexCentersX(k))^2 ) >4
                vortexCentersX_all(imNo,k) = NaN;
                vortexCentersY_all(imNo,k) = NaN;
                rotation_all(imNo,k) = NaN;
                vortexSize_all(imNo,k) = NaN;
                continue; % Move to the next vortex without breaking the loop
            end

            vortexSize_all(imNo,k) = sum(isolatedMask(:)) * pixel_area;  % Calculate the vortex size
            boundaries = bwboundaries(isolatedMask);
            boundaryIdx = boundaries{1};
            distances = nan(length(boundaryIdx), 1);
            vtan = nan(length(boundaryIdx), 1);
            theta = nan(length(boundaryIdx), 1);
    
            for idx = 1:length(boundaryIdx)
                % Get the row and column indices of the boundary pixel
                rIdx = boundaryIdx(idx, 1);
                cIdx = boundaryIdx(idx, 2);
                
                % Get physical coordinates for pixel (assumed from calibrated meshgrid)
                xp = x(rIdx, cIdx);
                yp = y(rIdx, cIdx);
                
                % Compute radial vector components relative to vortex center
                dx_val = xp - vortexCentersX(1, k);
                dy_val = yp - vortexCentersY(1, k);
                theta(idx) = atan2(dy_val, dx_val);    
                r = sqrt(dx_val^2 + dy_val^2);
                distances(idx) = r;
                
                if r == 0
                    continue;
                end
            
                % Get the velocity components at the boundary pixel
                up = u(rIdx, cIdx);
                vp = v(rIdx, cIdx);
                
                % Tangential velocity component
                vtan(idx) = -up * sin(theta(idx)) + vp * cos(theta(idx));
            end
            
            % Compute the angular velocity (assuming distances are in meters and converting to rad/s)
            rotation = vtan ./ (distances * 10^-3); % [rad/s]
            rotation_all(imNo,k) = mean(rotation / (2 * pi), 'omitnan');

        end
        % Identify opening columns for first_row_index
        openCols = find(~b_mask(first_row_index, :));
        if ~isempty(openCols)
            c_min = min(openCols);
            c_max = max(openCols);
        else
            c_min = NaN; 
            c_max = NaN;
        end
        
        % Define rows to evaluate
        rowOffsets = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5];        
        for idx = 1:numel(rowOffsets)
            r = first_row_index + rowOffsets(idx);  % Compute the actual row index
            
            if r >= 1 && r <= size(v,1)  % Ensure row is within bounds
                if rowOffsets(idx) < 0
                    % Use the same c_min, c_max as the reference row
                    if ~isnan(c_min)
                        v_flux_all(imNo,idx) = trapz(x(r, c_min:c_max), v(r, c_min:c_max));  % Compute v-flux
                        convective_velocity(imNo,idx) = mean(u(r, c_min:c_max), 'omitnan');  % Compute average u-velocity
                    end
                else
                    % Recalculate opening for each row below
                    newOpenCols = find(~b_mask(r, :));
                    if ~isempty(newOpenCols)
                        c_minBelow = min(newOpenCols);
                        c_maxBelow = max(newOpenCols);
                        v_flux_all(imNo,idx) = trapz(x(r, c_minBelow:c_maxBelow), v(r, c_minBelow:c_maxBelow));
                        convective_velocity(imNo,idx) = mean(u(r, c_minBelow:c_maxBelow), 'omitnan');  % Compute average u-velocity
                    end
                end
            end
        end


        if mod(imNo, 1000) == 0 || imNo ==1
            [~,dy] = gradient(Co_ords.Co_ords(5).y); %hard coded
            dy = -dy;
            dx       = gradient(Co_ords.Co_ords(5).x);
    
            swirl = swirling(u, v, dx, dy);           
    
            gammaField = gamma1(x, y, u, v, 10);
            gammaField(b_mask) = 0;
            figure(1);
            subplot(2,2,1);
            contourf(x, y, u, 20); colorbar;
            title(['u-velocity at t = ', num2str(imNo)]); xlabel('x'); ylabel('y');
            
            subplot(2,2,2);
            contourf(x, y, v, 20); colorbar; caxis([-0.1 0.1])
            title(['v-velocity at t = ', num2str(imNo)]); xlabel('x'); ylabel('y');
    
            subplot(2,2,3)
            imagesc(xcorners,ycorners,regionMask)
            title(['masked regions', num2str(imNo)]); xlabel('x'); ylabel('y');
            set(gca, 'YDir', 'normal')
           
            
            subplot(2,2,4);
            contourf(x, y, gammaField, 20, 'LineColor', 'none'); colorbar; caxis([-1 1]); hold on;
            plot(vortexCentersX_all(imNo,:), vortexCentersY_all(imNo,:), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
            title(['gamma1 Field and Detected Vortex Centres at t = ', num2str(imNo)]);
            xlabel('x'); ylabel('y');
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% 
            saveas(figure(1), fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated', sprintf('vortex_break_down_%05d.fig', imNo)));
            saveas(figure(1), fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated', sprintf('vortex_break_down_%05d.jpg', imNo)));
%             close(figure(1));
            close(gcf)
        end

        
    end

    segment_length = ceil(samplingRate(i) / 0.3);
    segment_length = min(caseImages, segment_length);
    overlap = segment_length / 2;
    nfft = 2^nextpow2(segment_length);
    frequencies = struct();
    
    % Define save directory
    save_dir = fullfile(base{i}, 'Statistics', num2str(imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    
    save_path = fullfile(save_dir, 'cavity_stats.mat');
    


    % Process X & Y Locations
    for k = 1:size(vortexCentersY_all, 2)
        frequencies.vortexCenterY(k) = process_and_plot(vortexCentersY_all(:, k), ['vortexCenterY_' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
        frequencies.vortexCenterX(k) = process_and_plot(vortexCentersX_all(:, k), ['vortexCenterX_' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
        frequencies.rotation(k) = process_and_plot(rotation_all(:,k), ['rotationRate' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
        frequencies.vortexSize(k) = process_and_plot(vortexSize_all(:,k), ['VortexSize' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
    
    end
    
    for k = 1:size(x_loc,2)
        frequencies.deltaStar(k) = process_and_plot(deltaStar_all(:,k), ['delta_star' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
        frequencies.theta(k) = process_and_plot(momentumThickness_all(:,k), ['Momentum thickness' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
    end
    
    
    for k = 1:size(v_flux_all,2)
        frequencies.vflux(k) = process_and_plot(v_flux_all(:,k), ['vflux' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
        frequencies.convective_velocity(k) = process_and_plot(convective_velocity(:,k), ['convective velocity' num2str(k)], save_dir, samplingRate(i), segment_length, overlap, nfft);
    
    
    end

        save(save_path, 'vortexCentersX_all', 'vortexCentersY_all', 'rotation_all', ...
    'v_flux_all', 'deltaStar_all', 'vortexSize_all', 'frequencies', 'momentumThickness_all', 'convective_velocity', 'POD_rebuild');


end

