function Inst_statistics(setup,Type,endpoint, use_merged_data)

    % /******************************************************************************
%  * Function: Perform_PIV_statistics_inst
%  * ---------------------------------------------------------------------------
%  * Description:
%  * This function calculates and saves statistical data (mean, RMS, fluctuation, 
%  * vorticity, divergence) from instantaneous PIV (Particle Image Velocimetry) 
%  * results. The function processes PIV data from a sequence of images for a 
%  * given camera and computes the following statistics:
%  *   - Mean velocity components (U, V)
%  *   - RMS velocity components (U, V)
%  *   - Velocity fluctuations (U', V')
%  *   - Vorticity
%  *   - Divergence of the velocity field
%  * The function operates for each camera in the setup and handles multiple runs 
%  * of PIV data.
%  * 
%  * For each camera, the function loads the PIV results, processes the data 
%  * frame by frame, and computes statistics such as the root mean square (RMS) 
%  * of velocity, the mean velocities, and the fluctuation components. It also 
%  * calculates the vorticity and divergence fields, and saves the statistical 
%  * data to a specified directory.
%  *
%  * The function uses parallel processing to speed up the computation by using 
%  * the `parfor` loop, distributing the image processing across multiple workers.
%  *
%  * ---------------------------------------------------------------------------
%  * Inputs:
%  *   struct setup : A structure containing all necessary configuration data
%  *                  and parameters for the PIV analysis, including:
%  *       - setup.imProperties.cameraCount  : The number of cameras.
%  *       - setup.imProperties.imageCount   : The number of images to process.
%  *       - setup.instantaneous.runs        : Number of runs for instantaneous 
%  *                                         PIV analysis.
%  *       - setup.instantaneous.nameConvention : Cell array defining the naming 
%  *                                             convention for the PIV data files.
%  *       - setup.instantaneous.windowSize  : The window size for each run.
%  *       - setup.directory.base            : The base directory for loading and 
%  *                                         saving data.
%  *       - setup.directory.statistics      : The base directory for saving statistics.
%  *   
%  *   Example:
%  *   setup = struct('imProperties', struct('cameraCount', 2, 'imageCount', 100),
%  *                  'instantaneous', struct('runs', 1:10, 'nameConvention', {...},
%  *                                         'windowSize', [32, 32]),
%  *                  'directory', struct('base', '/path/to/data', 'statistics', '/path/to/save/stats'));
%  *   
%  * ---------------------------------------------------------------------------
%  * Outputs:
%  *   - The function does not return values directly but saves the statistical 
%  *     data to files within the directory specified in the `setup` structure.
%  *     The output files are saved as `.mat` files in the directory:
%  *     `setup.directory.base/Statistics/<imageCount>/Cam<cameraNo>/Instantaneous/`
%  *     containing the following variables:
%  *       - mean_U         : Mean U velocity field.
%  *       - mean_V         : Mean V velocity field.
%  *       - RMS_U          : RMS of U velocity field.
%  *       - RMS_V          : RMS of V velocity field.
%  *       - U_prime_Uprime_mean  : Mean of U' * U' (velocity fluctuations).
%  *       - V_prime_Vprime_mean  : Mean of V' * V' (velocity fluctuations).
%  *       - U_prime_Vprime_mean  : Mean of U' * V' (cross-fluctuation).
%  *       - mean_Vorticity : Mean vorticity field.
%  *       - divergenceField: Divergence field.
%  *   
%  * ---------------------------------------------------------------------------
%  * Procedures:
%  *   1. For each camera, loop through each run specified in `setup.instantaneous.runs`.
%  *   2. For each run, load the PIV data and the coordinates from the respective 
%  *      directory.
%  *   3. Compute the mean and RMS velocity fields for both U and V.
%  *   4. For each image, calculate the velocity fluctuation fields (U' and V') 
%  *      and compute the cross-fluctuations and the second moments.
%  *   5. Calculate the vorticity and divergence fields for the velocity data.
%  *   6. Save the computed statistics in `.mat` files in the appropriate directory.
%  *   7. The `parfor` loop is used for parallel processing to speed up the computation.
%  *   
%  * ---------------------------------------------------------------------------
%  * Example Usage:
%  *   setup = struct('imProperties', struct('cameraCount', 1, 'imageCount', 50),
%  *                  'instantaneous', struct('runs', 1:5, 
%  *                                         'nameConvention', {'piv_data_%04d.mat'}, 
%  *                                         'windowSize', [32, 32]),
%  *                  'directory', struct('base', '/path/to/data', 'statistics', '/path/to/save/stats'));
%  *   
%  *   Perform_PIV_statistics_inst(setup);
%  *   
%  *   This will compute and save the statistical results for the PIV data across 
%  *   the specified cameras and runs, saving the results in the appropriate directory.
%  *   
%  *****************************************************************************/
    if setup.pipeline.statistics_inst
        fprintf('Inst_statistics for base: %s at %s\n', setup.directory.base, datetime('now'));
        
        % Loop through all cameras
        for CameraNo = 1:setup.imProperties.cameraCount
            fprintf('Processing Camera %d of %d\n', CameraNo, setup.imProperties.cameraCount);
            if use_merged_data
                dataloc = fullfile(setup.directory.base, 'Merged', type, endpoint);
            
            else      
                dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
               
                
            end
        
        for i = setup.instantaneous.runs
            VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
            Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
            b_mask = VelData.piv_result(i).b_mask;
            ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
            xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            dx       = gradient(Co_ords.Co_ords(i).x);
            [~,dy]   = gradient(Co_ords.Co_ords(i).y);
            dy = -dy;

            
            U_gridTotalsquares=zeros(size(VelData.piv_result(i).ux));
            V_gridTotalsquares=zeros(size(VelData.piv_result(i).uy));
            U_gridTotal = zeros(size(VelData.piv_result(i).ux));
            V_gridTotal = zeros(size(VelData.piv_result(i).uy));
            U_prime_Uprime_total = zeros(size(VelData.piv_result(i).ux));
            V_prime_Vprime_total = zeros(size(VelData.piv_result(i).uy));
            U_prime_Vprime_total = zeros(size(VelData.piv_result(i).ux));
            vorticity_total = zeros(size(VelData.piv_result(i).ux));
            divergenceField = zeros(size(VelData.piv_result(i).ux));
            parfor ImNo =1:setup.imProperties.imageCount
                VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, ImNo))]));
                U_gridTotalsquares = U_gridTotalsquares + (VelData.piv_result(i).ux).^2;
                V_gridTotalsquares = V_gridTotalsquares + (VelData.piv_result(i).uy).^2;
                U_gridTotal        = U_gridTotal + (VelData.piv_result(i).ux);
                V_gridTotal        = V_gridTotal + (VelData.piv_result(i).uy);
                [dvx,~]=gradient(VelData.piv_result(i).uy);
                [~,duy]=gradient(VelData.piv_result(i).ux);
                duy = -duy;
                vorticity_total = vorticity_total + (dvx./dx - duy./(dy));
                divergenceField= divergenceField + divergence(VelData.piv_result(i).ux, VelData.piv_result(i).uy);

            end
            RMS_U = sqrt(U_gridTotalsquares / setup.imProperties.imageCount);
            RMS_V = sqrt(V_gridTotalsquares / setup.imProperties.imageCount);
            mean_U = U_gridTotal / setup.imProperties.imageCount;
            mean_V = V_gridTotal / setup.imProperties.imageCount;
            divergenceField = divergenceField/setup.imProperties.imageCount;
                
            Setup_parpool(setup, 'Processes')
            parfor ImNo = 1:setup.imProperties.imageCount
                VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, ImNo))]));
                Ufluc=(VelData.piv_result(i).ux)-mean_U;
                Vfluc=(VelData.piv_result(i).uy)-mean_V;
                U_prime_Uprime_total=U_prime_Uprime_total+(Ufluc.*Ufluc);    
                V_prime_Vprime_total=V_prime_Vprime_total+(Vfluc.*Vfluc);
                U_prime_Vprime_total=U_prime_Vprime_total+(Vfluc.*Ufluc);
            end
            U_prime_Uprime_mean = U_prime_Uprime_total / setup.imProperties.imageCount;
            V_prime_Vprime_mean = V_prime_Vprime_total / setup.imProperties.imageCount;
            U_prime_Vprime_mean=U_prime_Vprime_total/setup.imProperties.imageCount;
            mean_Vorticity = vorticity_total/setup.imProperties.imageCount;
        
            if use_merged_data
                directory = fullfile(setup.directory.base, 'Statistics', 'Merged', Type, endpoint, 'meanStats');
            else
               
                directory=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], Type, endpoint, 'meanStats');
            end
            if ~exist(directory,'dir')
                mkdir(directory)
            end
            save(fullfile(directory,['MeanStats',num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)),'.mat']),"mean_U","mean_V","RMS_U","RMS_V",'U_prime_Uprime_mean', 'V_prime_Vprime_mean','U_prime_Vprime_mean', "mean_Vorticity","divergenceField")
        
            
            variableName ='Mean U ';
            plot_save_mask(mean_U, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
            variableName ='Mean V ';
            plot_save_mask(mean_V, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
            variableName ='RMS U ';
            plot_save_mask(RMS_U, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
            variableName ='RMS V ';
            plot_save_mask(RMS_V, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
    

            figure('Visible', 'off')
            fig = gcf;
            ax1 = axes;
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

            % Improved streamlines with better styling
            h = streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, mean_U, mean_V, 7);
            set(h, 'Color', 'k');
            title([Type, ' Streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2))], 'FontSize',setup.figures.titleFontSize,'Interpreter', 'latex');
            daspect([1 1 1]);
            fig = gcf;
            xlim([Co_ords.Co_ords(i).x(1,1) Co_ords.Co_ords(i).x(1,end)]);
            ylim([Co_ords.Co_ords(i).y(end,1) Co_ords.Co_ords(i).y(1,1)]);
            set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

            % Hold on to add mask layer
            hold(ax1, 'on');

            % Create mask overlay with light grey color (same as plot_save_mask)
            mask_vis = double(b_mask);
            [rows, cols] = size(mask_vis);
            mask_rgb = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8); % Light grey RGB

            % Display the mask as an RGB image overlay
            j = image(ax1, xcorners, ycorners, mask_rgb);
            set(j, 'AlphaData', mask_vis * 0.7); % Show masked areas with transparency

            daspect(ax1, [1 1 1]);
            hold(ax1, 'off');

            saveas(fig, fullfile(directory, [ '_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg']));
            saveas(fig, fullfile(directory, ['_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc']));
            saveas(fig, fullfile(directory, ['_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig']));
            close(fig);
            close all;
            delete(j);


            [ny, nx] = size(mean_U);
            if ny <= nx
                square_size = ny;
                x_start = floor((nx - ny)/2) + 1;
                x_end = x_start + square_size - 1;
                y_start = 1;
                y_end = ny;
            else
                square_size = nx;
                y_start = floor((ny - nx)/2) + 1;
                y_end = y_start + square_size - 1;
                x_start = 1;
                x_end = nx;
            end

            % Crop all arrays to the central square
            U_cropped = mean_U(y_start:y_end, x_start:x_end);
            V_cropped = mean_V(y_start:y_end, x_start:x_end);
            b_mask_cropped = b_mask(y_start:y_end, x_start:x_end);
            vorticity_cropped = mean_Vorticity(y_start:y_end, x_start:x_end);
            x_co_ords_cropped = Co_ords.Co_ords(i).x(y_start:y_end, x_start:x_end);
            y_co_ords_cropped = Co_ords.Co_ords(i).y(y_start:y_end, x_start:x_end);
            xcorners_cropped = [x_co_ords_cropped(1,1), x_co_ords_cropped(1,end)];
            ycorners_cropped = [y_co_ords_cropped(1,1), y_co_ords_cropped(end, end)];

            % Upsample
            upsample_factor = 5;
            upsampled_mean_u = imresize(U_cropped, upsample_factor, "bicubic");
            upsampled_mean_v = imresize(V_cropped, upsample_factor, "bicubic");
            upsampled_b_mask = imresize(b_mask_cropped, upsample_factor, "nearest");
            upsampled_vorticity = imresize(vorticity_cropped, upsample_factor, "bicubic");

            upsample_cat_v = cat(3, -upsampled_mean_v, upsampled_mean_u);
            upsample_cat_v = perform_vf_normalization(upsample_cat_v);

            options.spot_size = 2;
            options.flow_correction = 1;

            lic = perform_lic(upsample_cat_v, 12, options);
            lic(upsampled_b_mask==1) = 0;
            upsampled_vorticity(upsampled_b_mask==1) = 0;

            % Mask LIC and vorticity arrays
            lic_masked = lic;
            lic_masked(upsampled_b_mask==1) = NaN;


            % Plot with improved style (like plot_save_mask)
            figure('Visible', 'off');
            fig = gcf;
            set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            ax1 = axes;
            hold(ax1, 'on');

            % 1. LIC background (masked)
            A = imagesc(ax1, xcorners_cropped, ycorners_cropped, lic_masked);
            colormap(ax1, gray);
            set(A, 'AlphaData', ~isnan(lic_masked));

            % 2. Mask overlay (draw last, flip vertically for correct orientation with YDir normal)
            mask_vis = double(upsampled_b_mask);
            [rows, cols] = size(mask_vis);
            mask_rgb = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
            j = image(ax1, xcorners_cropped, ycorners_cropped, (mask_rgb));
            set(j, 'AlphaData', (mask_vis) * 0.7);

            % Axis and colorbar settings
            daspect(ax1, [1 1 1]);
            set(ax1, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, ...
                'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
            title( 'Line Integral Convolutions', ...
                'FontSize', setup.figures.titleFontSize, 'Interpreter', 'latex');

            hold(ax1, 'off');
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

            saveas(gcf, fullfile(directory, [ '_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) ...
                'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg']));
            saveas(gcf, fullfile(directory, [ '_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) ...
                'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc']));
            saveas(gcf, fullfile(directory, ['_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) ...
                'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig']));
            close(gcf)





            variableName ="u prime squared ";
            plot_save_mask(U_prime_Uprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="v prime squared ";
            plot_save_mask(V_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="u prime v prime ";
            plot_save_mask(U_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');            variableName ="Mean Vorticity ";
            plot_save_mask(mean_Vorticity, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="Mean Divergence ";
            plot_save_mask(U_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
        
            close all;
            

        end % End of runs loop
        end % End of camera loop
    end % End of pipeline check
end % End of function