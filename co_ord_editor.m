function co_ord_editor(setup, CameraNo)
    dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Instantaneous'));
    
    if setup.pipeline.instantaneous_cords    % Process instantaneous runs
        for i = setup.instantaneous.runs
            VelData = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
            Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
            
            % Get current coordinates
            ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
            xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            
            % Create interactive figure
            fig = figure('Name', ['Run ' num2str(i) ' - Use toolbar to zoom/pan, then click to set origin']);
            
            % Create masked data for visualization
            data_to_show = VelData.piv_result(i).ux;
            nan_mask = VelData.piv_result(i).b_mask;
            
            % Create grey background for NaN regions (same as plot_save_mask)
            [rows, cols] = size(nan_mask);
            grey_bg = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
            
            % Display data first
            h = imagesc(xcorners, ycorners, data_to_show);
            set(h, 'AlphaData', ~nan_mask);
            set(gca, 'YDir', 'normal');
            hold on;
            
            % Add grey overlay for NaN regions
            h_grey = image(xcorners, ycorners, grey_bg);
            set(h_grey, 'AlphaData', double(nan_mask) * 0.7);
            
            colorbar;
            title(['Run ' num2str(i) ' - Use zoom/pan tools, then press ENTER and click to set origin']);
            xlabel('X coordinate');
            ylabel('Y coordinate');
            axis equal;
            hold off;
            
            % Wait for user interaction with better instructions
            fprintf('\n=== Run %d Coordinate Editor ===\n', i);
            fprintf('1. Use the zoom and pan tools in the figure toolbar to navigate\n');
            fprintf('2. When ready to set origin, press ENTER in the command window\n');
            fprintf('3. Then click once on the plot to set the new origin (0,0)\n\n');
            
            % Wait for user to press enter before enabling click
            input('Press ENTER when ready to set origin, then click on the plot: ');
            
            % Now get the click
            [click_x, click_y] = ginput(1);
            
            % Close the interactive figure
            close(fig);
            
            % Calculate offset from clicked point
            x_offset = click_x;
            y_offset = click_y;
            
            % Update coordinates by subtracting the offset
            Co_ords.Co_ords(i).x = Co_ords.Co_ords(i).x - x_offset;
            Co_ords.Co_ords(i).y = Co_ords.Co_ords(i).y - y_offset;
            
            % Save updated coordinates
            save(fullfile(dataloc, 'Co_ords.mat'), '-struct', 'Co_ords');
            
            % Create confirmation plot with new coordinates
            fig_confirm = figure('Name', ['Run ' num2str(i) ' - Updated coordinates (non-interactive)']);
            new_ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
            new_xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            
            % Display updated data
            h_confirm = imagesc(new_xcorners, new_ycorners, data_to_show);
            set(h_confirm, 'AlphaData', ~nan_mask);
            set(gca, 'YDir', 'normal');
            hold on;
            
            % Add grey overlay for NaN regions in confirmation
            h_grey_confirm = image(new_xcorners, new_ycorners, grey_bg);
            set(h_grey_confirm, 'AlphaData', double(nan_mask) * 0.7);
            
            colorbar;
            title(['Run ' num2str(i) ' - Coordinates updated! Origin at clicked point']);
            xlabel('X coordinate (updated)');
            ylabel('Y coordinate (updated)');
            axis equal;
            
            % Add crosshairs at origin
            plot(0, 0, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
            plot([-max(abs(new_xcorners)), max(abs(new_xcorners))], [0, 0], 'r--', 'LineWidth', 1);
            plot([0, 0], [-max(abs(new_ycorners)), max(abs(new_ycorners))], 'r--', 'LineWidth', 1);
            hold off;
            
            % Display for 3 seconds then close
            pause(3);
            close(fig_confirm);
            
            fprintf('Coordinates updated for run %d. New origin set at clicked location.\n', i);
            
            % Wait for user acknowledgment before proceeding to next run
            if i < max(setup.instantaneous.runs)
                input('Press Enter to continue to next run...');
            end
        end
        
        fprintf('Coordinate editing complete for all instantaneous runs.\n');
    end 
    % Process ensemble runs
    if setup.pipeline.ensemble_cords
       
        dataloc_ensemble = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Ensemble'));
        
        for i = setup.ensemble.runs
            VelData = load(fullfile(dataloc_ensemble, [num2str(sprintf(setup.ensemble.nameConvention{1}, 1))]));
            Co_ords = load(fullfile(dataloc_ensemble, 'Co_ords.mat'));
            
            % Get current coordinates
            ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
            xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            
            % Create interactive figure
            fig = figure('Name', ['Ensemble Run ' num2str(i) ' - Use toolbar to zoom/pan, then click to set origin']);
            
            % Create masked data for visualization
            data_to_show = VelData.piv_result(i).ux;
            nan_mask = isnan(data_to_show);
            
            % Create grey background for NaN regions
            [rows, cols] = size(nan_mask);
            grey_bg = cat(3, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8, ones(rows, cols) * 0.8);
            
            % Display data first
            h = imagesc(xcorners, ycorners, data_to_show);
            set(h, 'AlphaData', ~nan_mask);
            set(gca, 'YDir', 'normal');
            hold on;
            
            % Add grey overlay for NaN regions
            h_grey = image(xcorners, ycorners, grey_bg);
            set(h_grey, 'AlphaData', double(nan_mask) * 0.7);
            
            colorbar;
            title(['Ensemble Run ' num2str(i) ' - Use zoom/pan tools, then press ENTER and click to set origin']);
            xlabel('X coordinate');
            ylabel('Y coordinate');
            axis equal;
            hold off;
            
            % Wait for user interaction with better instructions
            fprintf('\n=== Ensemble Run %d Coordinate Editor ===\n', i);
            fprintf('1. Use the zoom and pan tools in the figure toolbar to navigate\n');
            fprintf('2. When ready to set origin, press ENTER in the command window\n');
            fprintf('3. Then click once on the plot to set the new origin (0,0)\n\n');
            
            % Wait for user to press enter before enabling click
            input('Press ENTER when ready to set origin, then click on the plot: ');
            
            % Now get the click
            [click_x, click_y] = ginput(1);
            
            % Close the interactive figure
            close(fig);
            
            % Calculate offset from clicked point
            x_offset = click_x;
            y_offset = click_y;
            
            % Update coordinates by subtracting the offset
            Co_ords.Co_ords(i).x = Co_ords.Co_ords(i).x - x_offset;
            Co_ords.Co_ords(i).y = Co_ords.Co_ords(i).y - y_offset;
            
            % Save updated coordinates
            save(fullfile(dataloc_ensemble, 'Co_ords.mat'), '-struct', 'Co_ords');
            
            % Create confirmation plot with new coordinates
            fig_confirm = figure('Name', ['Ensemble Run ' num2str(i) ' - Updated coordinates (non-interactive)']);
            new_ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
            new_xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            
            % Display updated data
            h_confirm = imagesc(new_xcorners, new_ycorners, data_to_show);
            set(h_confirm, 'AlphaData', ~nan_mask);
            set(gca, 'YDir', 'normal');
            hold on;
            
            % Add grey overlay for NaN regions in confirmation
            h_grey_confirm = image(new_xcorners, new_ycorners, grey_bg);
            set(h_grey_confirm, 'AlphaData', double(nan_mask) * 0.7);
            
            colorbar;
            title(['Ensemble Run ' num2str(i) ' - Coordinates updated! Origin at clicked point']);
            xlabel('X coordinate (updated)');
            ylabel('Y coordinate (updated)');
            axis equal;
            
            % Add crosshairs at origin
            plot(0, 0, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
            plot([-max(abs(new_xcorners)), max(abs(new_xcorners))], [0, 0], 'r--', 'LineWidth', 1);
            plot([0, 0], [-max(abs(new_ycorners)), max(abs(new_ycorners))], 'r--', 'LineWidth', 1);
            hold off;
            
            % Display for 3 seconds then close
            pause(3);
            close(fig_confirm);
            
            fprintf('Coordinates updated for ensemble run %d. New origin set at clicked location.\n', i);
            
            % Wait for user acknowledgment before proceeding to next run
            if i < max(setup.ensemble.runs)
                input('Press Enter to continue to next ensemble run...');
            end
        end
        
        fprintf('Coordinate editing complete for all runs.\n');
    end
end