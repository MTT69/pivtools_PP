function POST_POD_UV_Batches(setup, Type, CameraNo, endpoint, domain)

    if setup.pipeline.POST_POD_batches

        if strcmp(Type, 'Calibrated')        
            dataloc = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ...
                               ['Cam', num2str(CameraNo)], 'Instantaneous', endpoint);
        elseif strcmp(Type, 'Uncalibrated')
            dataloc = fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ...
                               ['Cam', num2str(CameraNo)], 'Instantaneous');
        else
            error('Incompatible type for instantaneous statistics');
        end

        fprintf('Performing POD analysis using SVD at %s\n', datetime('now'));

        for i = setup.instantaneous.runs
            VelData = load(fullfile(dataloc, sprintf(setup.instantaneous.nameConvention{1}, 1)));
            Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
            ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
            xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
            
            mask = VelData.piv_result(i).b_mask;
            row_average = mean(mask, 2);
            index = find(row_average > 0.8, 1);

            if strcmp(domain, 'Above')
                mask(index:end, :) = 1; 
            elseif strcmp(domain, 'Below')
                se = strel('disk', 3);  % Structuring element for dilation
                mask = imdilate(mask, se);  % Create buffer zone
                mask(1:index-1, :) = 1;

            end

            h = size(Co_ords.Co_ords(i).y, 1);
            w = size(Co_ords.Co_ords(i).x, 2);
            for batch = 1 : (setup.imProperties.imageCount/setup.imProperties.caseImages)

            
                % Preallocate U and V matrices
                U = nan( h * w, setup.imProperties.caseImages);
                V = nan( h * w, setup.imProperties.caseImages);

                
    
                % Loop over the number of images per batch.
                parfor idx = 1 : setup.imProperties.caseImages
                    % Compute the actual image number from the batch and the loop index.
                    ImNo = (batch-1)*setup.imProperties.caseImages + idx;
                    
                    % Load the data for each image.
                    VelData = load(fullfile(dataloc, sprintf(setup.instantaneous.nameConvention{1}, ImNo)));
                    u = VelData.piv_result(i).ux;
                    v = VelData.piv_result(i).uy;
                
                    % Apply the mask.
                    u(mask) = 0;
                    v(mask) = 0;
                
                    % Store the velocity fields using the loop index as the sliced index.
                    U(:, idx) = u(:);
                    V(:, idx) = v(:);
                end

                
                fprintf('POD Vectors Loaded at %s\n', datetime('now'));
                mean_u = mean(U, 2);
                mean_v = mean(V,2);
    
                % Subtract mean flow from U and V
                U = U - mean_u;
                V = V - mean_v;
    
                % Concatenate U and V vertically
                UV = [U; V];  % Size: (2*Nt x Nx)
    
                % Apply SVD (economy mode)
                [U_svd, S_svd, V_svd] = svd(UV, 'econ');
    
                % Extract singular values and compute energy contribution
                singular_values = diag(S_svd);
                totalEnergy = sum(singular_values.^2);
                energyPortion = (singular_values.^2 / totalEnergy) * 100;
                cumulativeEnergy = cumsum(energyPortion);
    
                % Find number of modes capturing 80% of energy
                targetPercentage = 80;
                numModes = find(cumulativeEnergy >= targetPercentage, 1);
                numModes = min(numModes, 10); % Limit to 10 modes
    
                % Normalize spatial modes
                PHI = normc(V_svd);
    
                % Save results
                directory = fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ...
                                     ['Cam' num2str(CameraNo)], 'Instantaneous', Type, endpoint, domain,num2str(batch));
                if ~exist(directory, 'dir')
                    mkdir(directory);
                end
    
                filename = sprintf('POD_stats_%dx%d.mat', w, h);
                save(fullfile(directory, filename), 'PHI', 'U_svd', 'S_svd', 'V_svd', ...
                     'singular_values', 'energyPortion', 'cumulativeEnergy', 'h', 'w','mean_u', 'mean_v', "index", '-v7.3');
    
                % Plot cumulative energy
                figure;
                plot(1:min(setup.imProperties.imageCount, 50), cumulativeEnergy(1:min(setup.imProperties.imageCount, 50)), 'b-o', 'LineWidth', 2, 'MarkerSize', 10);
                ylabel('Cumulative Energy (%)', 'FontSize', setup.figures.labelFontSize);
                xlabel('Mode Index', 'FontSize', setup.figures.labelFontSize);
                title(['Pass ', num2str(i), ' Cumulative Energy Distribution'], 'FontSize', setup.figures.titleFontSize, 'Interpreter', 'latex');
                grid on; box on;
                set(gca, 'FontSize', setup.figures.axisFontSize);
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
                % Save figure
                saveas(gcf, fullfile(directory, ['CumulativeEnergyDistribution_pass', num2str(i), '.jpg']));
                saveas(gcf, fullfile(directory, ['CumulativeEnergyDistribution_pass', num2str(i), '.epsc']));
                saveas(gcf, fullfile(directory, ['CumulativeEnergyDistribution_pass', num2str(i), '.fig']));
                close(gcf);
    
    

                for k = 1:numModes
                    Mode_spatial_U = reshape(U_svd(1:h*w, k), [h, w]); % First h*w rows for U spatial mode
                    Mode_spatial_V = reshape(U_svd(h*w+1:end, k), [h, w]); % Last h*w rows for V spatial mode
                    Mode_spatial_U(mask > 0) = 0;
                    Mode_spatial_V(mask > 0) = 0;
                    variableName = ['Pass ', num2str(i), ' U Mode ', num2str(k),' ' ];
                    plot_save_mask(Mode_spatial_U, mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], i, variableName, 'New','dir')
                    openfig(fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig'],"")));
                    hold on
                    streamPlot = streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, Mode_spatial_U, Mode_spatial_V, 7);
                    set(streamPlot, 'Color', 'k'); % Set the color of the streamslice to black
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg'],"")))
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc'],"")))
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig'],"")))
                    close(gcf)
                end
    
                % Plot V modes
                for k = 1:numModes
                    Mode_spatial_U = reshape(U_svd(1:h*w, k), [h, w]); % First h*w rows for U spatial mode
                    Mode_spatial_V = reshape(U_svd(h*w+1:end, k), [h, w]); % Last h*w rows for V spatial mode
                    variableName = ['Pass ', num2str(i), ' V Mode ', num2str(k),' ' ];
                    plot_save_mask(Mode_spatial_V, mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,directory,setup.instantaneous.runs,[setup.instantaneous.windowSize(i,1), setup.instantaneous.windowSize(i,2)], i, variableName, 'New','dir')
                    openfig(fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig'],"")));
                    hold on
                    streamPlot = streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, Mode_spatial_U, Mode_spatial_V, 7);
                    set(streamPlot, 'Color', 'k'); % Set the color of the streamslice to black
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.jpg'],"")))
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.epsc'],"")))
                    saveas(gcf, fullfile(directory,join([variableName, num2str(setup.instantaneous.windowSize(i,1)) 'x' num2str(setup.instantaneous.windowSize(i,2)) '.fig'],"")))
                    close(gcf)
                end
            end
        end
    end
end
