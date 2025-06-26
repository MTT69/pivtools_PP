function SumStatistics(setup,Type,endpoint, use_merged_data)
% SumStatistics calculates and saves sum statistics for PIV data.
% Type: must be 'Ensemble'
% endpoint: string for merged/endpoint folder
% use_merged_data: boolean, if true, process merged data (Type must match merged folder)
%
% The function processes PIV data and saves sum statistics for each run.
% Output directories and file loading are handled based on use_merged_data.

    if nargin < 4
        use_merged_data = false;
    end

    if setup.pipeline.statistics_sum
        % Only allow 'Ensemble' type
        if ~strcmp(Type, 'Ensemble')
            error('SumStatistics only supports Type = ''Ensemble''.');
        end

        for CameraNo = 1:setup.imProperties.cameraCount
            if use_merged_data
                dataloc = fullfile(setup.directory.base, 'Merged', Type, endpoint);
                folderPath = fullfile(setup.directory.base, 'Statistics', 'Merged', Type, endpoint);
            else
                dataloc = fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount), ['Cam', num2str(CameraNo)], 'Ensemble',endpoint);
                folderPath = fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Ensemble', endpoint);
            end

            if ~exist(folderPath, 'dir')
                mkdir(folderPath)
            end

            % Main processing loop
            for i = setup.ensemble.runs
                data = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
                b_mask = data.piv_result(i).b_mask;
                Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));

                ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)];
                xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];

                variableName = "u prime squared ";
                plot_save_mask(data.piv_result(i).Uturb, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'Sum', instfolderPath);

                variableName = "v prime squared ";
                plot_save_mask(data.piv_result(i).Vturb, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'Sum', instfolderPath);

                variableName = "u prime v prime ";
                plot_save_mask(data.piv_result(i).UturbVturb, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'Sum', instfolderPath);

                variableName = "Mean U ";
                plot_save_mask(data.piv_result(i).ux, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'Sum', instfolderPath);

                variableName = "Mean V ";
                plot_save_mask(data.piv_result(i).uy, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'Sum', instfolderPath);


                % % Optionally, add Uncalibrated/MergeUncalibrated support if needed
                % if strcmp(Type, 'Uncalibrated') || strcmp(Type, 'MergeUncalibrated')

                %     variableName = "Peak Heights AB";
                %     plot_save_mask(data.piv_result(i).peakheights_AB, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sxy Cross Correlation";
                %     plot_save_mask(data.piv_result(i).sig_AB_xy, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sx Cross Correlation";
                %     plot_save_mask(data.piv_result(i).sig_AB_x, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sy Cross Correlation";
                %     plot_save_mask(data.piv_result(i).sig_AB_y, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sxy Auto Correlation";
                %     plot_save_mask(data.piv_result(i).sig_A_xy, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sx Auto Correlation";
                %     plot_save_mask(data.piv_result(i).sig_A_x, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                %     variableName = "Sy Auto Correlation";
                %     plot_save_mask(data.piv_result(i).sig_A_y, b_mask, xcorners, ycorners, setup.figures.axisFontSize, setup.figures.titleFontSize, folderPath, setup.ensemble.runs, [setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)], i, variableName, 'New', 'dir');

                % end
            end
        end
    end
end