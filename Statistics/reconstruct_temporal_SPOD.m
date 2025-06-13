function reconstruct_temporal_SPOD(setup,CameraNo,endpoint,combine, modes, frequencies,nDFT,run)

    if setup.pipeline.SPOD

        directory = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Instantaneous',endpoint,'SPOD');
        SPOD_stats = load(fullfile(directory, 'SPOD_results.mat'));
        
        % f = (0:nDFT-1)/dt/nDFT;
        % reducing A reduces the mode set
        A_reduced = SPOD_stats.A;
        A_reduced(:, modes, :) = 0;  
        A_reduced(frequencies, :, :) = 0;
        nOvlp = nDFT / 2;
        reconstructedFrames = invspod(SPOD_stats.P, A_reduced, nDFT, nOvlp);
        temporalMean_full = repmat(SPOD_stats.temporalMean, [1, 1, setup.imProperties.imageCount]);
        temporalMean_full = permute(temporalMean_full, [3, 1, 2]);
        save_directory = fullfile(directory,['Reconstruction modes ' num2str(modes(1)), ' _ ' num2str(modes(end)), ' Frequencies ' num2str(SPOD_stats.f(frequencies(1))), ' _ ', num2str(SPOD_stats.f(frequencies(2))) 'hz']);
        if ~exist(save_directory, 'dir')
            mkdir(save_directory);
        end
        if combine
            if height(reconstructedFrames) ~= height(temporalMean_full)
                disp(' Not all frames can be reconstructed as DFT is not a multiple of nImages')
            end
            for imNo = 1: height(reconstructedFrames)
                piv_result = struct();
                
                file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
                file_path = fullfile(save_directory, file_name);
                piv_result(run).ux = squeeze(reconstructedFrames(imNo,1:end/2,:))+ squeeze(temporalMean_full(imNo,1:end/2,:));
                piv_result(run).uv = squeeze(reconstructedFrames(imNo,end/2+1:end,:)) + squeeze(temporalMean_full(imNo, end/2+1:end,:));
                PIVSAVE(file_path,piv_result)
            end
        else
            if height(reconstructedFrames) ~= height(temporalMean_full)
                disp(' Not all frames can be reconstructed as DFT is not a multiple of nImages')
            end
            for imNo = 1: height(reconstructedFrames)
                piv_result = struct();
                file_name = num2str(sprintf(setup.instantaneous.nameConvention{1}, imNo));
                file_path = fullfile(save_directory, file_name);
                piv_result(run).ux = squeeze(reconstructedFrames(imNo,:,:))+ squeeze(temporalMean_full(imNo,:,:));
                PIVSAVE(file_path,piv_result)
            end
        end
    end
end