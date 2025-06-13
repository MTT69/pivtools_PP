function reconstructSPOD_modes(setup,CameraNo,endpoint,combine, modes, frequencies)
    if setup.pipeline.SPOD
        directory = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous',endpoint,'SPOD');
        SPOD_stats = load(fullfile(directory, 'SPOD_results.mat'));
        for fi = frequencies
            for mi = modes
                if combine
                    % First plot (U velocity)
                    figure(1);
                    imagesc(real(squeeze(SPOD_stats.P(fi, 1:end/2, :, mi)))), axis equal tight;
                    title(['U velocity f=' num2str(SPOD_stats.f(fi), '%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(SPOD_stats.L(fi, mi), '%.2g')]);
                    
                    % Save the figure
                    saveas(gcf, fullfile(directory, ['U_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.jpg']));
                    saveas(gcf, fullfile(directory, ['U_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.epsc']));
                    saveas(gcf, fullfile(directory, ['U_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.fig']));
                    close(gcf); % Close the figure
                
                    % Second plot (V velocity)
                    figure(1);
                    imagesc(real(squeeze(SPOD_stats.P(fi, end/2+1:end, :, mi)))), axis equal tight;
                    title(['V velocity f=' num2str(SPOD_stats.f(fi), '%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(SPOD_stats.L(fi, mi), '%.2g')]);
                    
                    % Save the figure
                    saveas(gcf, fullfile(directory, ['V_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.jpg']));
                    saveas(gcf, fullfile(directory, ['V_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.epsc']));
                    saveas(gcf, fullfile(directory, ['V_velocity_f_' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.fig']));
                    close(gcf); % Close the figure
                else
                    % Combined plot
                    figure(1);
                    imagesc(real(squeeze(SPOD_stats.P(fi, :, :, mi)))), axis equal tight;
                    title(['f=' num2str(SPOD_stats.f(fi), '%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(SPOD_stats.L(fi, mi), '%.2g')]);
                    
                    % Save the figure
                    saveas(gcf, fullfile(directory, ['U' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.jpg']));
                    saveas(gcf, fullfile(directory, ['U' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.epsc']));
                    saveas(gcf, fullfile(directory, ['U' num2str(SPOD_stats.f(fi), '%.2f') '_mode_' num2str(mi) '.fig']));
                    close(gcf); % Close the figure
                end

        
            end
        end
    end
end