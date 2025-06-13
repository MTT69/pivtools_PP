function Inst_statistics(setup,Type,CameraNo,endpoint)

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
        if strcmp(Type,'Calibrated')        
            dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous',endpoint));
            
        elseif strcmp(Type, 'Uncalibrated')
            dataloc = (fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
        else
            error('Incompatible type for instantaneous statistics');
            
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
            for k = 1:numel(mean_U)
                if mean_U(k)<0
                    RMS_U(k) = -RMS_U(k);
                end
                if mean_V(k) <0
                    RMS_V(k) = - RMS_V(k);
                end
            end
                
            Setup_parpool(setup, 'Processes')
            parfor ImNo = 1:setup.imProperties.imageCount
                VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, ImNo))]));
                Ufluc=(VelData.piv_result(i).ux)-RMS_U;
                Vfluc=(VelData.piv_result(i).uy)-RMS_V;
                U_prime_Uprime_total=U_prime_Uprime_total+(Ufluc.*Ufluc);    
                V_prime_Vprime_total=V_prime_Vprime_total+(Vfluc.*Vfluc);
                U_prime_Vprime_total=U_prime_Vprime_total+(Vfluc.*Ufluc);
            end
            U_prime_Uprime_mean = U_prime_Uprime_total / setup.imProperties.imageCount;
            V_prime_Vprime_mean = V_prime_Vprime_total / setup.imProperties.imageCount;
            U_prime_Vprime_mean=U_prime_Vprime_total/setup.imProperties.imageCount;
            mean_Vorticity = vorticity_total/setup.imProperties.imageCount;
        

            directory=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', Type, endpoint);
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
            streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y,(mean_U),(mean_V), 7)
            title([Type, ' Streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2))], 'FontSize',setup.figures.titleFontSize,'Interpreter', 'latex');
            daspect([1 1 1]);
            fig = gcf;
            xlim([Co_ords.Co_ords(i).x(1,1) Co_ords.Co_ords(i).x(1,end)]);
            ylim([Co_ords.Co_ords(i).y(end,1) Co_ords.Co_ords(i).y(1,1)]);
            set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels
            mask = double(b_mask);
            mask(mask==1) =inf;
            ax2 = axes;
            j = imagesc(xcorners,ycorners, mask);
            set(j, 'AlphaData', isinf(mask));
            set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels            colormap(ax2, [0 0 0]);
            daspect(ax2, [1 1 1]);
            
            % Align axes
            linkaxes([ax1, ax2]);
            
            % Hide the second set of axes
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            
            % Set proper positioning of axes
            ax2.Position = ax1.Position;
            
            % Make sure the figure aspect ratio is consistent
            daspect([1 1 1]);
            
            % Remove extra space between plots
            ax1.XLim = ax2.XLim;
            ax1.YLim = ax2.YLim;
            
            saveas(fig, fullfile(directory, [ '_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg']));
            saveas(fig, fullfile(directory, ['_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc']));
            saveas(fig, fullfile(directory, ['_Mean_streamlines ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig']));
            close(fig);
            close all;
            delete(j);


%             %TODO only works for square imagery
%             figure('Visible', 'off')
%             center_col = floor(size(mean_U, 2) / 2);
%             half_width = floor(size(mean_U, 1) / 2);
%             U_cropped = mean_U(:, center_col - half_width : center_col + half_width);
%             V_cropped = mean_V(:, center_col - half_width : center_col + half_width);
%             x_co_ords_cropped = Co_ords.Co_ords(i).x(:, center_col - half_width : center_col + half_width);
%             y_co_ords_cropped = Co_ords.Co_ords(i).y(:, center_col - half_width : center_col + half_width);
%             xcorners_cropped = [x_co_ords_cropped(1,1), x_co_ords_cropped(1,end)];
%             ycorners_cropped = [y_co_ords_cropped(1,1), y_co_ords_cropped(end, end)];
%             b_mask_cropped = b_mask(:, center_col - half_width : center_col + half_width);
%             vorticity_cropped = mean_Vorticity(:, center_col - half_width : center_col + half_width);
%             upsampled_mean_u = imresize(U_cropped, 5, "bicubic");
%             upsampled_mean_v = imresize(V_cropped, 5, "bicubic");
%             upsampled_b_mask = imresize(b_mask_cropped, 5, "nearest");
%             upsampled_vorticity = imresize(vorticity_cropped, 5, "bicubic");
% 
%             upsample_cat_v = cat(3, -upsampled_mean_v, upsampled_mean_u);
%             upsample_cat_v = perform_vf_normalization(upsample_cat_v);
% 
% 
% 
%             options.spot_size = 2;
%             options.flow_correction = 1;
% 
%             lic = perform_lic(upsample_cat_v, 12, options);
%             lic(upsampled_b_mask==1) =0;
%             upsampled_vorticity(upsampled_b_mask==1)=0;
% 
%             ax1 = axes;
%             A = imagesc(xcorners_cropped, ycorners_cropped,lic);
%             colormap(ax1, gray);
%             daspect(ax1, [1 1 1]);
%             alphaData_lic = lic ~= 0; % True for non-zero values, false for zeros
%             set(A, 'AlphaData', alphaData_lic);
%             set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels            colormap(ax2, [0 0 0]);
%             hold on;
%             ax2 = axes;
%             lower_limit = prctile(upsampled_vorticity(:),5);
%             upper_limit = prctile(upsampled_vorticity(:),95);
%             B = imagesc(xcorners_cropped, ycorners_cropped,upsampled_vorticity,[lower_limit, upper_limit]);
%             colormap(ax2, redbluezero(lower_limit, upper_limit));
%             daspect(ax2, [1 1 1]);
%             colorbar;
%             alphaData = ~isnan(upsampled_vorticity);
%             set(B, 'AlphaData', alphaData * 0.5); % Use 0.5 transparency for non-zero values
%             set(gca, 'YDir', 'normal', 'FontSize', setup.figures.axisFontSize, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex'); % Thicker axes, LaTeX labels            colormap(ax2, [0 0 0]);
% 
%             % % Align axes
%             linkaxes([ax1, ax2]);
%             % 
%             % % Hide the second set of axes
%             ax2.Visible = 'off';
%             ax2.XTick = [];
%             ax2.YTick = [];
%             ax2.HitTest = 'off'; % Ignore mouse clicks on ax2
%             ax2.PickableParts = 'none'; % Prevent ax2 from capturing interactions
% 
%             ax2.Position = ax1.Position;
%             daspect([1 1 1]);
%             % 
%             % % Remove extra space between plots
% 
%             ax1.XLim = ax2.XLim;
%             ax1.YLim = ax2.YLim;
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% 
%             title([Type, '_Vorticity_LIC ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2))], 'FontSize',setup.figures.titleFontSize,'Interpreter', 'latex');
%             saveas(gcf, fullfile(directory, [ '_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg']));
%             saveas(gcf, fullfile(directory, [ '_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc']));
%             saveas(gcf, fullfile(directory, ['_vorticity_lic ' num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig']));
% 
%             close(gcf)





            variableName ="u prime squared ";
            plot_save_mask(U_prime_Uprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="v prime squared ";
            plot_save_mask(V_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="u prime v prime ";
            plot_save_mask(U_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="Mean Vorticity ";
            plot_save_mask(mean_Vorticity, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');

            variableName ="Mean Divergence ";
            plot_save_mask(U_prime_Vprime_mean, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)],i, variableName, 'Inst','dir');
        
            close all;
            

        end
    end
end