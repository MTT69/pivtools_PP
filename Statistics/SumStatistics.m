function SumStatistics(setup,Type,CameraNo)

    if setup.pipeline.statistics_sum
        if strcmp(Type,'Calibrated')        
            dataloc = (fullfile(setup.directory.base, 'CalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Ensemble'));
            
        elseif strcmp(Type, 'Uncalibrated')
            dataloc = (fullfile(setup.directory.base, 'UncalibratedPIV', num2str(setup.imProperties.imageCount),['Cam', num2str(CameraNo)], 'Ensemble'));
        else
            error('Incompatible type for Ensemble statistics');
            
        end

        folderPath=fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'Ensemble', Type);
        if ~exist(folderPath,'dir')
            mkdir(folderPath)
        end

        instfolderPath =fullfile(setup.directory.base, 'Statistics', num2str(setup.imProperties.imageCount), ['Cam' num2str(CameraNo)], 'instantaneous', Type);

        for i = setup.ensemble.runs
            data = load(fullfile(dataloc, [num2str(sprintf(setup.instantaneous.nameConvention{1}, 1))]));
            b_mask = data.piv_result(i).b_mask;
            Co_ords = load(fullfile(dataloc,'Co_ords.mat'));
            

            ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; 
            xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];

            variableName ="u prime squared " ;
            plot_save_mask(data.piv_result(i).Uturb, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'Sum',instfolderPath);
            
            variableName ="v prime squared ";
            plot_save_mask(data.piv_result(i).Vturb, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'Sum', instfolderPath);
            
            variableName ="u prime v prime ";
            plot_save_mask(data.piv_result(i).UturbVturb, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'Sum', instfolderPath);

            variableName ="Mean U ";
            plot_save_mask(data.piv_result(i).ux, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'Sum', instfolderPath);
            
            variableName ="Mean V ";
            plot_save_mask(data.piv_result(i).uy, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'Sum', instfolderPath);
            
        
            if strcmp(Type, 'Uncalibrated') || strcmp(Type, 'MergeUncalibrated')

                variableName ="Peak Heights AB";
                plot_save_mask(data.piv_result(i).peakheights_AB, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');
            
                variableName ="Sxy Cross Correlation";
                plot_save_mask(data.piv_result(i).sig_AB_xy, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');
                
                variableName ="Sx Cross Correlation";
                plot_save_mask(data.piv_result(i).sig_AB_x, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');

                variableName ="Sy Cross Correlation";
                plot_save_mask(data.piv_result(i).sig_AB_y, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');

                variableName ="Sxy Auto Correlation";
                plot_save_mask(data.piv_result(i).sig_A_xy, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');
                
                variableName ="Sx Auto Correlation";
                plot_save_mask(data.piv_result(i).sig_A_x, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');

                variableName ="Sy Auto Correlation";
                plot_save_mask(data.piv_result(i).sig_A_y, b_mask,xcorners,ycorners,setup.figures.axisFontSize,setup.figures.titleFontSize,folderPath,setup.ensemble.runs,[setup.ensemble.windowSize(i,1), setup.ensemble.windowSize(i,2)],i, variableName, 'New', 'dir');
                
            end
        end
    end
end