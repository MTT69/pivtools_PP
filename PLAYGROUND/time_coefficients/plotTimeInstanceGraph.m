  function plotTimeInstanceGraph( b_mask, xcorners_top, ycorners_top, field_top, lower_top, upper_top, ...
        xcorners_bottom, ycorners_bottom, field_bottom, lower_bottom, upper_bottom, ...
        time_coefficient, imNo, dt, xlim_l, xlim_u, ylim_l, ylim_u, plot_type,index)
    

    
        % Create a new figure
        fig = figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
        % Mask creation
        mask = double(b_mask);
        mask(mask == 1) = inf;
        ax3 = axes;
        j = imagesc(xcorners_bottom,ycorners_bottom, mask(index:end,:));
        set(j, 'AlphaData', isinf(mask(index:end,:)));
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

        % center_col = floor(size(ux, 2) / 2);
        % half_width = floor(size(ux,1) / 2);
        % U_cropped = ux(:, center_col - half_width : center_col + half_width);
        % V_cropped = uy(:, center_col - half_width : center_col + half_width);
        % b_mask_bottom = b_mask(index:end,:);
        % b_mask_cropped = b_mask_bottom(:, center_col - half_width : center_col + half_width);
        % scale_factor = 5;
        % U_upsampled = imresize(U_cropped, scale_factor, 'bicubic');
        % V_upsampled = imresize(V_cropped, scale_factor, 'bicubic');
        % b_mask_upsampled = imresize(b_mask_cropped, scale_factor, 'bicubic');
        % U_upsampled(b_mask_upsampled==1) =0;
        % V_upsampled(b_mask_upsampled==1) =0;
        % v_upsampled = cat(3, -V_upsampled, U_upsampled);
        % v_upsampled = perform_vf_normalization(v_upsampled);
        % options.spot_size = 2;
        % options.flow_correction = 1;

        % lic_result_upsampled = perform_lic(v_upsampled, 12, options);
        % lic_result_upsampled(b_mask_upsampled==1) =0;

        % ax4 = axes;
        % img = imagesc(ax4, xcorners_bot_cropped, ycorners_bot_cropped, lic_result_upsampled);
        % colormap(ax4, 'gray');
        % daspect(ax4, [1 1 1]);
        % set(img, 'AlphaData', 0.4); % Set transparency to 50%
        % set(ax4, 'YDir', 'normal', 'FontSize', 18);
        % ax4.Position = ax2.Position;
        % linkaxes([ax4, ax3]);
        % ax4.XLim = ax2.XLim;
        % ax4.YLim = ax2.YLim;
        % ax4.XTick = [];
        % ax4.YTick = [];
        % ax4.Visible = 'off';
    
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
        else
            new_position = [0.05, 0.05, 0.9, 0.9]; % Nearly full-screen position
            ax1.Position = new_position;
            ax2.Position = new_position;
            ax3.Position = new_position;

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