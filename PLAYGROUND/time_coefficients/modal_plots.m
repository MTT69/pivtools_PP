
base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
    'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
    'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};
imageCount = 18000;
caseImages = 3600;
CameraNo =1;
domain = 'Below';
i =5;
windowSize=  [128 128; 64 64; 32 32; 16 16; 16 16];
runs =[5];
numModes =10;


for dir = 1:length(base)
    for batch = 1 : (imageCount/caseImages)

        directory = fullfile(base{dir}, 'Statistics', num2str(imageCount), ...
            ['Cam' num2str(CameraNo)], 'Instantaneous', 'Calibrated', domain,num2str(batch));
        dataloc = fullfile(base{dir}, 'CalibratedPIV', num2str(imageCount), ...
            ['Cam', num2str(CameraNo)], 'Instantaneous');

        VelData = load(fullfile(dataloc, [num2str(sprintf('%05d.mat', 1))]));
        Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
        ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
        xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];

        mask = VelData.piv_result(i).b_mask;
        row_average = mean(mask, 2);
        index = find(row_average > 0.8, 1);

        if strcmp(domain, 'Above')
            mask(index:end, :) = 1; 
        elseif strcmp(domain, 'Below')
            mask(1:index-1, :) = 1;
            se = strel('disk', 3);  % Structuring element for dilation
            mask = imdilate(mask, se);  % Create buffer zone
        end

        load(fullfile(directory,'POD_stats_319x149.mat'))

        for k = 1:numModes
            Mode_spatial_U = reshape(U_svd(1:h*w, k), [h, w]); % First h*w rows for U spatial mode
            Mode_spatial_V = reshape(U_svd(h*w+1:end, k), [h, w]); % Last h*w rows for V spatial mode
            Mode_spatial_U(mask > 0) = 0;
            Mode_spatial_V(mask > 0) = 0;
            variableName = ['Pass ', num2str(i), ' U Mode ', num2str(k),' ' ];
            plot_save_mask(Mode_spatial_U, mask,xcorners,ycorners,16,20,directory,runs,[windowSize(i,1), windowSize(i,2)], i, variableName, 'New','dir')
            openfig(fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.fig'],"")));
            hold on
            streamPlot = streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, Mode_spatial_U, Mode_spatial_V, 7);
            set(streamPlot, 'Color', 'k'); % Set the color of the streamslice to black
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.jpg'],"")))
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.epsc'],"")))
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.fig'],"")))
            close(gcf)
        end
        
        % Plot V modes
        for k = 1:numModes
            Mode_spatial_U = reshape(U_svd(1:h*w, k), [h, w]); % First h*w rows for U spatial mode
            Mode_spatial_V = reshape(U_svd(h*w+1:end, k), [h, w]); % Last h*w rows for V spatial mode
            variableName = ['Pass ', num2str(i), ' V Mode ', num2str(k),' ' ];
            plot_save_mask(Mode_spatial_V, mask,xcorners,ycorners,16,20,directory,runs,[windowSize(i,1), windowSize(i,2)], i, variableName, 'New','dir')
            openfig(fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.fig'],"")));
            hold on
            streamPlot = streamslice(Co_ords.Co_ords(i).x, Co_ords.Co_ords(i).y, Mode_spatial_U, Mode_spatial_V, 7);
            set(streamPlot, 'Color', 'k'); % Set the color of the streamslice to black
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.jpg'],"")))
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.epsc'],"")))
            saveas(gcf, fullfile(directory,join([variableName, num2str(windowSize(i,1)) 'x' num2str(windowSize(i,2)) '.fig'],"")))
            close(gcf)
        end
    end
end