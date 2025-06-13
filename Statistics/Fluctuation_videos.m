function Fluctuation_videos(setup, endpoint, u_inf, sampling_rate, plot_type, domain, modes, dividing_row, frame_range, CameraNo, Type)


    if setup.pipeline.Fluctuation_videos
        disp(['Running fluctuating videos for base ', setup.directory.base]);
        for i = setup.instantaneous.runs
            for k = modes
                directory_path = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', endpoint);
                file_name = fullfile(directory_path,num2str(sprintf(setup.instantaneous.nameConvention{1}, 1)));
                piv = load(file_name);
                ux = piv.piv_result(i).ux; % X-direction velocity
                [h,w] = size(ux);
                
                dt = 1/sampling_rate;
                lower_top = -0.05; upper_top = 0.05;
                lower_bottom = -0.025; upper_bottom = 0.025;
                stats_directory=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', Type, endpoint, domain);
                filename = ['POD_stats_', num2str(w), 'x', num2str(h), '.mat'];
                Data = load(fullfile(stats_directory, filename));
                A = Data. A;
                clear Data
                frames = frame_range*2+1;
                xlim_l = -80;
                xlim_u = 80;
                ylim_u = 30;
                ylim_l = -45;
                for p_t_loop = 1:2
                    if p_t_loop == 1
                        p_t = 'peak';
                    else
                        p_t = 'trough';
                    end
                    mp4_filename = fullfile(stats_directory, ['FlowDevelopment-', p_t, '-', plot_type, '-Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mp4']);

                    if strcmp(p_t, 'peak')
                        peak_data = load(fullfile(stats_directory, ['peak_indices_Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']));
                        combined_indices = peak_data.extended_peak_locations;
                        indices = peak_data.peak_locations;
                    elseif strcmp(p_t, 'trough')
                        trough_data = load(fullfile(stats_directory, ['trough_indices_Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']));
                        combined_indices = trough_data.extended_trough_locations;
                        indices = trough_data.trough_locations;
                    end
                    if strcmp(plot_type, 'u')
                        mean_data =load(fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', endpoint, ['MeanStats', num2str(setup.instantaneous.windowSize(i,1)), 'x', num2str(setup.instantaneous.windowSize(i,2)), '.mat']));
                        mean_field = mean_data.mean_U/ u_inf;
                    elseif strcmp(plot_type, 'v')
                        mean_data =load(fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', endpoint,['MeanStats', num2str(setup.instantaneous.windowSize(i,1)), 'x', num2str(setup.instantaneous.windowSize(i,2)), '.mat']));
                        mean_field = mean_data.mean_V/ u_inf;
                    end
                    Co_ords = load(fullfile(directory_path, "Co_ords"));
                    ycorners_top = [Co_ords.Co_ords(i).y(1, 1), Co_ords.Co_ords(i).y(dividing_row, end)];
                    xcorners_top = [Co_ords.Co_ords(i).x(1, 1), Co_ords.Co_ords(i).x(dividing_row, end)];
                    center_col = floor(size(ux(dividing_row:end,:), 2) / 2);
                    half_width = floor(size(ux(dividing_row:end,:),1) / 2);
                    Co_ords_bottom_cropped_y = Co_ords.Co_ords(i).y(dividing_row:end,center_col - half_width : center_col + half_width);
                    Co_ords_bottom_cropped_x = Co_ords.Co_ords(i).x(dividing_row:end, center_col - half_width : center_col + half_width);
                    % Upsample the coordinates
                    scale_factor = 5;
                    Co_ords_bottom_cropped_y_upsampled = imresize(Co_ords_bottom_cropped_y, scale_factor, 'bicubic');
                    Co_ords_bottom_cropped_x_upsampled = imresize(Co_ords_bottom_cropped_x, scale_factor, 'bicubic');
                    ycorners_bot_cropped = [Co_ords_bottom_cropped_y_upsampled(1, 1), Co_ords_bottom_cropped_y_upsampled(end, end)];                    
                    xcorners_bot_cropped = [Co_ords_bottom_cropped_x_upsampled(1, 1), Co_ords_bottom_cropped_x_upsampled(end, end)];

                    ycorners_bottom = [Co_ords.Co_ords(i).y(dividing_row, 1), Co_ords.Co_ords(i).y(end, end)];
                    xcorners_bottom = [Co_ords.Co_ords(i).x(dividing_row, 1), Co_ords.Co_ords(i).x(end, end)];
                    % used for scalar so full resolution
                    Total_flow_top = zeros(size(ux(1:dividing_row,:)));
                    Total_flow_bottom = zeros(size(ux(dividing_row:end,:)));
                    Total_ux_top = zeros(size(ux(1:dividing_row,:)));
                    Total_uy_top = zeros(size(ux(1:dividing_row,:)));
                    Total_ux_bottom = zeros(size(ux(dividing_row:end,:)));
                    Total_uy_bottom = zeros(size(ux(dividing_row:end,:)));

                    numCells = length(indices);
                    % This is preallocation for the shorter peak average video
                    developmentVideo_top = cell(1, numCells);
                    developmentVideo_bottom = cell(1, numCells);
                    developmentVideo_ux_top = cell(1, numCells);
                    developmentVideo_ux_bottom = cell(1, numCells);
                    developmentVideo_uy_top = cell(1, numCells);
                    developmentVideo_uy_bottom = cell(1, numCells);
                    % Initialize each cell with the relevant size array

                    for H = 1:numCells
                        developmentVideo_top{H} = zeros(size(Total_flow_top)); % Match Total_flow_top dimensions
                        developmentVideo_bottom{H} = zeros(size(Total_flow_bottom)); % Match Total_flow_bottom dimensions
                        developmentVideo_ux_top{H} = zeros(size(Total_flow_top)); % Match Total_flow_top dimensions
                        developmentVideo_uy_top{H} = zeros(size(Total_flow_top)); % Match Total_flow_top dimensions
                        developmentVideo_ux_bottom{H} = zeros(size(Total_flow_bottom)); % Match Total_flow_bottom dimensions
                        developmentVideo_uy_bottom{H} = zeros(size(Total_flow_bottom)); % Match Total_flow_bottom dimensions
                        
                    end
                    v = VideoWriter(mp4_filename, 'MPEG-4');  % Create a VideoWriter object for MP4
                    v.FrameRate = 3;  % Set the frame rate (e.g., 10 frames per second)
                    open(v);  % Open the video writer
                    count = 0;
                    index = 0;
                    for imNo = combined_indices
                        index = index+1;
                        file_name = fullfile(directory_path,num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo)));
                        piv=load(file_name);
                        ux = piv.piv_result(i).ux/ u_inf;
                        uy = piv.piv_result(i).uy/ u_inf;
                        if strcmp(plot_type, 'u')
                            fluc = ux - mean_field;
                        elseif strcmp(plot_type, 'v')
                            fluc = uy - mean_field;
                        end
                        b_mask = piv.piv_result(i).b_mask;
                        % total fluctating
                        fluc_top = fluc(1:dividing_row,:);
                        fluc_bottom = fluc(dividing_row:end,:);
                        ux_top = ux(1:dividing_row,:);
                        uy_top = uy(1:dividing_row,:);
                        ux_bottom = ux(dividing_row:end,:);
                        uy_bottom = uy(dividing_row:end,:);
                        Total_flow_top = Total_flow_top + fluc_top;
                        Total_flow_bottom = Total_flow_bottom + fluc_bottom;
                        Total_ux_top = Total_ux_top + ux_top;
                        Total_uy_top = Total_uy_top + uy_top;
                        Total_ux_bottom = Total_ux_bottom + ux_bottom;
                        Total_uy_bottom = Total_uy_bottom + uy_bottom;
                        % Average accross peaks for secondary video
                        if mod(index - (frames+1) ,frames) == 0 && imNo >= frames
                            count = count +1;
                        end
                        developmentVideo_bottom{count} = developmentVideo_bottom{count} + fluc_bottom;
                        developmentVideo_top{count} = developmentVideo_top{count} + fluc_top;
                        developmentVideo_ux_top{count} = developmentVideo_ux_top{count} + ux_top;
                        developmentVideo_uy_top{count} = developmentVideo_uy_top{count} + uy_top;
                        developmentVideo_ux_bottom{count} = developmentVideo_ux_bottom{count} + ux_bottom;
                        developmentVideo_uy_bottom{count} = developmentVideo_uy_bottom{count} + uy_bottom;
                        
                        time_coefficient = A(1:setup.imProperties.imageCount, k);

                        plotTimeInstanceGraph(xcorners_bot_cropped, ycorners_bot_cropped, ux_bottom,uy_bottom,b_mask, xcorners_top, ycorners_top, fluc_top, lower_top, upper_top, ...
                        xcorners_bottom, ycorners_bottom, fluc_bottom, lower_bottom, upper_bottom,...
                        time_coefficient, imNo, dt, xlim_l, xlim_u, ylim_l, ylim_u, 'Fluctuating',dividing_row);

                        % Capture the current frame
                        frame = getframe(gcf);

                        % Write the frame to the video
                        writeVideo(v, frame);

                        % Close the current figure
                        close(gcf);
                    end
                    close(v);
                    % average fluctuating component
                    Average_flow_bottom = Total_flow_bottom/length(combined_indices); 
                    Average_flow_top = Total_flow_top/length(combined_indices);
                    % average velocity component
                    Average_ux_bottom = Total_ux_bottom/length(combined_indices);
                    Average_uy_bottom = Total_uy_bottom/length(combined_indices);
                    Average_ux_top = Total_ux_top/length(combined_indices);
                    Average_uy_top = Total_uy_top/length(combined_indices);
                  
    
    
                   

                    plotTimeInstanceGraph(xcorners_bot_cropped, ycorners_bot_cropped,Average_ux_bottom,Average_uy_bottom,...
                    piv.piv_result(i).b_mask, xcorners_top, ycorners_top, Average_flow_top, lower_top, upper_top, ...
                    xcorners_bottom, ycorners_bottom, Average_flow_bottom, lower_bottom, upper_bottom, ...
                    [], [], [], xlim_l, xlim_u, ylim_l, ylim_u, 'average',dividing_row);

                    file_name = fullfile(stats_directory, ['AverageFlow-', p_t, '-', plot_type, '-Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mat']);
                    % Save the figure as a JPEG file
                    saveas(gcf, [file_name, '.jpeg']);

                    % Save the figure as a MATLAB figure file
                    savefig(gcf, [file_name, '.fig']);
                    close gcf
                    mp4_filename = fullfile(stats_directory, ['FlowDevelopment-average', p_t, '-', plot_type, '-Mode', num2str(k), '_', num2str(w), 'x', num2str(h), '.mp4']);
                    

                    v = VideoWriter(mp4_filename, 'MPEG-4');  % Create a VideoWriter object for MP4
                    v.FrameRate = 3;  % Set the frame rate (e.g., 10 frames per second)
                    open(v);  % Open the video writer
                    
                    for H = 1:numCells
                        developmentVideo_top{H} = developmentVideo_top{H} / frames;
                        developmentVideo_bottom{H} = developmentVideo_bottom{H} / frames;
                        developmentVideo_ux_top{H} = developmentVideo_ux_top{H} / frames;
                        developmentVideo_uy_top{H} = developmentVideo_uy_top{H} / frames;
                        developmentVideo_ux_bottom{H} = developmentVideo_ux_bottom{H} / frames;
                        developmentVideo_uy_bottom{H} = developmentVideo_uy_bottom{H} / frames;
                        plotTimeInstanceGraph(xcorners_bot_cropped, ycorners_bot_cropped,developmentVideo_ux_bottom{H}, developmentVideo_uy_bottom{H}, b_mask, xcorners_top, ycorners_top, developmentVideo_top{H}, lower_top, upper_top, ...
                            xcorners_bottom, ycorners_bottom, developmentVideo_bottom{H}, lower_bottom, upper_bottom, ...
                            time_coefficient, indices(H), dt, xlim_l, xlim_u, ylim_l, ylim_u,'Fluctuating',dividing_row);
                        
                        % Capture frame
                        frame = getframe(gcf);
                        writeVideo(v, frame);

                        % Close figure
                        close(gcf);
                    end
                    close(v)

                end

                

            end
        end
    end
    function plotTimeInstanceGraph(xcorners_bot_cropped, ycorners_bot_cropped,ux, uy, b_mask, xcorners_top, ycorners_top, field_top, lower_top, upper_top, ...
        xcorners_bottom, ycorners_bottom, field_bottom, lower_bottom, upper_bottom, ...
        time_coefficient, imNo, dt, xlim_l, xlim_u, ylim_l, ylim_u, plot_type,dividing_row)
    

    
        % Create a new figure
        fig = figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
        % Mask creation
        mask = double(b_mask);
        mask(mask == 1) = inf;
        ax3 = axes;
        j = imagesc(xcorners_bottom,ycorners_bottom, mask(dividing_row:end,:));
        set(j, 'AlphaData', isinf(mask(dividing_row:end,:)));
        set(gca, 'YDir', 'normal', 'FontSize', 18);
        colormap(ax3, [0 0 0]);
    
        % Align axes
        daspect([1 1 1]);
        ax3.XTick = [];
        ax3.YTick = [];
        ax3.Visible = 'off';
    
        % Plot bottom field data
        ax2 = axes;
        imagesc(ax2, xcorners_bottom, ycorners_bottom, field_bottom, [lower_bottom, upper_bottom]);
        daspect(ax2, [1 1 1]);
        set(ax2, 'YDir', 'normal', 'FontSize', 18);
        colormap(ax2, redbluezero(lower_bottom, upper_bottom));
        ax3.Position = ax2.Position;
        uistack(ax2, 'bottom'); % Push ax2 to the bottom of the stack
        linkaxes([ax2, ax3]);
        ax3.XLim = ax2.XLim;
        ax3.YLim = ax2.YLim;
        cbar2 = colorbar(ax2, 'eastoutside');
        cbar2.Label.String = [plot_type ]; % Add label to colorbar if needed

        center_col = floor(size(ux, 2) / 2);
        half_width = floor(size(ux,1) / 2);
        U_cropped = ux(:, center_col - half_width : center_col + half_width);
        V_cropped = uy(:, center_col - half_width : center_col + half_width);
        b_mask_bottom = b_mask(dividing_row:end,:);
        b_mask_cropped = b_mask_bottom(:, center_col - half_width : center_col + half_width);
        scale_factor = 5;
        U_upsampled = imresize(U_cropped, scale_factor, 'bicubic');
        V_upsampled = imresize(V_cropped, scale_factor, 'bicubic');
        b_mask_upsampled = imresize(b_mask_cropped, scale_factor, 'bicubic');
        U_upsampled(b_mask_upsampled==1) =0;
        V_upsampled(b_mask_upsampled==1) =0;
        v_upsampled = cat(3, -V_upsampled, U_upsampled);
        v_upsampled = perform_vf_normalization(v_upsampled);
        options.spot_size = 2;
        options.flow_correction = 1;

        lic_result_upsampled = perform_lic(v_upsampled, 12, options);
        lic_result_upsampled(b_mask_upsampled==1) =0;

        ax4 = axes;
        img = imagesc(ax4, xcorners_bot_cropped, ycorners_bot_cropped, lic_result_upsampled);
        colormap(ax4, 'gray');
        daspect(ax4, [1 1 1]);
        set(img, 'AlphaData', 0.4); % Set transparency to 50%
        set(ax4, 'YDir', 'normal', 'FontSize', 18);
        ax4.Position = ax2.Position;
        linkaxes([ax4, ax3]);
        ax4.XLim = ax2.XLim;
        ax4.YLim = ax2.YLim;
        ax4.XTick = [];
        ax4.YTick = [];
        ax4.Visible = 'off';
    
        % Plot top field data
        ax1 = axes;
        imagesc(ax1, xcorners_top, ycorners_top, field_top, [lower_top, upper_top]);
        daspect(ax1, [1 1 1]);
        set(ax1, 'YDir', 'normal', 'FontSize', 18);
        cbar1 = colorbar(ax1, 'westoutside');
        cbar1.Label.String = [plot_type ]; % Add label to colorbar if needed
        ax2.Position = ax1.Position;
        ax1.Visible = 'off';
        linkaxes([ax1, ax2, ax3]);



        if strcmp(plot_type, 'Fluctuating') && ~isempty(time_coefficient)
            new_position = [0.1, 0.4, 0.8, 0.55];
            ax1.Position = new_position;
            ax2.Position = new_position;
            ax3.Position = new_position;
            ax4.Position = new_position;
        else
            new_position = [0.05, 0.05, 0.9, 0.9]; % Nearly full-screen position
            ax1.Position = new_position;
            ax2.Position = new_position;
            ax3.Position = new_position;
            ax4.Position = new_position;

        end
        uistack(ax3, 'top'); % Push ax2 to the bottom of the stack
%         hold on;
%     
%         % Add quiver plots for velocity vectors
%         quiver(x_top, y_top, ux_top*2 , uy_top*2 , 'off', 'k'); % Top region
%         quiver(x_bottom, y_bottom, ux_bottom*15 , uy_bottom*15 , 'off', 'k'); % Bottom region
    
        xlim([xlim_l, xlim_u]);
        ylim([ylim_l, ylim_u]);
    
        % Add temporal coefficients plot if applicable
        if strcmp(plot_type, 'Fluctuating') && ~isempty(time_coefficient)
            ax5 = axes('Position', [0.1, 0.05, 0.8, 0.25]);
            plot((1:length(time_coefficient)) * dt, time_coefficient, 'LineWidth', 1.5);
            hold on;
    
            % Mark the current time instance
            current_time = imNo * dt;
            plot(current_time, time_coefficient(imNo), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
            hold off;
            xlabel('Time (s)');
            ylabel('Amplitude');
            title('Temporal Coefficients (Current Instance Marked)');
            grid on;
        end
    
        hold off;
    end
end