function  Correlation_stats(setup)
    
   
    if setup.pipeline.statistics_correlation
        fprintf('Performing Correlation_analysis at %s\n',(datetime('now')));
        for run = setup.instantaneous.runs
            all_CorMean = cell(1, setup.imProperties.cameraCount ); % averaged across instances
            all_MeanMat = cell(1, setup.imProperties.cameraCount ); % averaged in space
    
            DisplacementU = cell(1, setup.imProperties.cameraCount ); % averaged in space
            DisplacementV = cell(1, setup.imProperties.cameraCount ); % averaged in space
            NanCountAll   = cell(1,setup.imProperties.cameraCount );
            DataSize = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', sprintf( setup.instantaneous.nameConvention{1},1)));
            [x,y]=size(DataSize.piv_result(run).ux);
            for CameraNo=1:setup.imProperties.cameraCount 
                BufferMatCor = zeros(x, y);
                MeanMat = zeros(setup.imProperties.imageCount,1);
                MeanDispU = zeros(x,y);
                MeanDispV = zeros(x,y);
                NanCount = zeros(x,y);
                parfor ImNo = 1:setup.imProperties.imageCount
                    PIV = load(fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous', sprintf(setup.instantaneous.nameConvention{1},ImNo)));
                    pk_mag = PIV.piv_result(run).peak_mag;% PIV.piv_result(run)=  (finest) interrogation window results
                    dispU=PIV.piv_result(run).ux;
                    dispV=PIV.piv_result(run).uy;
                    pk_mag(PIV.piv_result(run).b_mask)=nan; % make masked out areas nans so they dont affect averages
                    dispU(PIV.piv_result(run).b_mask)=nan;
                    dispV(PIV.piv_result(run).b_mask)=nan;
                    NanCount = NanCount + PIV.piv_result(run).nan_mask;
                    MeanDispU=MeanDispU+dispU;
                    MeanDispV=MeanDispV+dispV;
                    BufferMatCor = BufferMatCor + pk_mag;% Taking good data
                    nonNanIndices = ~isnan(pk_mag);
                    MeanMat(ImNo) = sum(pk_mag(nonNanIndices), 'all') / sum(nonNanIndices, 'all'); % average across all good data in space
                end
                CorMean = BufferMatCor / (setup.imProperties.imageCount); % average across all snapshots i.e. time
                MeanDispU=MeanDispU/(setup.imProperties.imageCount);
                MeanDispV=MeanDispV/(setup.imProperties.imageCount);
                DisplacementU{CameraNo}=MeanDispU;
                DisplacementV{CameraNo}=MeanDispV;
                all_CorMean{CameraNo} = CorMean;
                all_MeanMat{CameraNo} = MeanMat;
                NanCount = 100*(NanCount/setup.imProperties.imageCount);
                NanCountAll{CameraNo} = NanCount;
                
                % Define the directory path
                folderPath = fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(1)], 'Instantaneous', 'Corr');
                
                % Ensure the folder exists
                if ~exist(folderPath, 'dir')
                    mkdir(folderPath);
                end
                
                % Save the file
                save(fullfile(folderPath, [num2str(run), '.mat']), "CorMean", "MeanMat", "NanCountAll");
    
                
            end
    
        
    
            saveLocation=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount));
            if ~exist(saveLocation, 'dir')
                mkdir(saveLocation);
            end
            
            for CameraNo = 1:setup.imProperties.cameraCount 
                subplot(3, setup.imProperties.cameraCount , CameraNo);
                imagesc(all_CorMean{CameraNo});
                title(['Camera ', num2str(CameraNo), ' CorMean']);
                colorbar
                clim([0, 1])
                pbaspect([1 1 1])
            end
            
            for CameraNo = 1:setup.imProperties.cameraCount 
                subplot(3, setup.imProperties.cameraCount , setup.imProperties.cameraCount  + CameraNo);
                plot(all_MeanMat{CameraNo})
                title(['Camera ', num2str(CameraNo), ' MeanMat']);
                ylim([0 1])
                pbaspect([1 1 1])
            end
    
            for CameraNo = 1: setup.imProperties.cameraCount 
                subplot(3, setup.imProperties.cameraCount , 2*setup.imProperties.cameraCount +CameraNo);
                imagesc(NanCountAll{CameraNo});
                clim([0,100])
                colorbar
                pbaspect([1 1 1])
                title(['Camera ', num2str(CameraNo), ' Temporal Nan percentage']);
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
            saveas(gcf, fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Correlation Evaluation', num2str(run), '.jpg']));
    
            saveas(gcf, fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Correlation Evaluation', num2str(run), '.fig']));
    
            close(gcf)
            close all;
    
            figure;
            % Initialize variables to store standard deviations for each camera
            P5_values = zeros(setup.imProperties.cameraCount , 1);
            P95_values = zeros(setup.imProperties.cameraCount , 1);
            
            % Loop through each camera
            for CameraNo = 1:setup.imProperties.cameraCount 
                % Extract the data for the current camera
                displacement_data = DisplacementU{CameraNo}(:);
                
                % Calculate the 5th percentile (P5) and 95th percentile (P95)
                P5_values(CameraNo) = prctile(displacement_data, 5);
                P95_values(CameraNo) = prctile(displacement_data, 95);
            end       
        
        
            for CameraNo =1:setup.imProperties.cameraCount 
                subplot(2,setup.imProperties.cameraCount , CameraNo);
                numBins = abs(floor(P5_values(CameraNo)-P95_values(CameraNo))*10); % this should give 1 bin per tenth of a m/s
                histogram(DisplacementU{CameraNo} , 'Normalization', 'pdf', 'BinWidth', 5, 'NumBins', numBins);
                xlabel('Pixel Displacement U');
                ylabel('Probability Density');
                title(['Camera ', num2str(CameraNo),' Peak Locking check'])
                pbaspect([1 1 1])
                xlim([P5_values(CameraNo), P95_values(CameraNo)]);
            end
            
            P5_values = zeros(setup.imProperties.cameraCount , 1);
            P95_values = zeros(setup.imProperties.cameraCount , 1);
            
            % Loop through each camera
            for CameraNo = 1:setup.imProperties.cameraCount 
                % Extract the data for the current camera
                displacement_data = DisplacementV{CameraNo}(:);
                
                % Calculate the 5th percentile (P5) and 95th percentile (P95)
                P5_values(CameraNo) = prctile(displacement_data, 5);
                P95_values(CameraNo) = prctile(displacement_data, 95);
            end
        
            for CameraNo =1:setup.imProperties.cameraCount 
                subplot(2,setup.imProperties.cameraCount , setup.imProperties.cameraCount +CameraNo);
                numBins = abs(floor(P5_values(CameraNo)-P95_values(CameraNo))*10); % this should give 1 bin per tenth of a m/s
                histogram(DisplacementV{CameraNo} , 'Normalization', 'pdf', 'BinWidth', 5, 'NumBins', numBins);
                xlabel('Pixel Displacement v');
                ylabel('Probability Density');
                title(['Camera ', num2str(CameraNo),' Peak Locking check'])
                ylim([0 1])
                pbaspect([1 1 1])
                xlim([P5_values(CameraNo), P95_values(CameraNo)]);
            
            end
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
            saveas(gcf, fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['PeakLocking Statistics', num2str(run), '.jpg']));
    
            saveas(gcf, fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['PeakLocking Statistics', num2str(run), '.fig']));
            close(gcf)
            close all;
        end
    end


end
      
     



