
clear variables;

CameraNo = 1;
base = { ...
    'D:\Full\Processed_PIV_validation\90degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\90degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\60degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\60degree_250light_250hz_1000dt', ...
%     'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt_reverse', ...
%     'D:\Full\Processed_PIV_validation\30degree_400light_100hz_3000dt', ...
%     'D:\Full\Processed_PIV_validation\30degree_250light_250hz_1000dt', ...
};
sampling_rate =[100,250,100,250,100,250,100,250,100,250,]; %hz 
u_inf = [0.34,1.18,0.34,1.18,0.34,1.18,0.34,1.18,0.34,1.18];
dt_piv = {3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6),3000*10^(-6),1000*10^(-6)}; % basic calibration dt - cell for different runs
imageCount = 18000;
caseImages = 3600;
i = 5;
plot_type ='v';
modes = [1:30];
scaleFactor = 9.53;
endpoint = "";

frequencies = [3,3];



for run = 1:length(base)
    dataloc = fullfile(base{run}, 'CalibratedPIV', num2str(imageCount), ...
    ['Cam', num2str(CameraNo)], 'Instantaneous', endpoint);
    statsloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
        ['Cam' num2str(CameraNo)], 'Instantaneous', ...
        'Calibrated');
    VelData = load(fullfile(dataloc, sprintf('%05d.mat', 1)));
    Co_ords = load(fullfile(dataloc, 'Co_ords.mat'));
    ycorners = [Co_ords.Co_ords(i).y(1,1), Co_ords.Co_ords(i).y(end, end)]; %, xcoords(m, 1), xcoords(m, n)];
    xcorners = [Co_ords.Co_ords(i).x(1,1), Co_ords.Co_ords(i).x(end, end)];
    [h, w] = size(VelData.piv_result(i).ux);
    SPOD = load(fullfile(statsloc, ['SPOD_stats_', num2str(w), 'x', num2str(h), '.mat']));

    A_reduced = zeros(size(SPOD.A));  

    A_reduced(frequencies(1):frequencies(2), modes(1):modes(2), :) = SPOD.A(frequencies(1):frequencies(2), modes(1):modes(2), :);

    % Reconstruct frames using only the selected frequency's modes
    reconstructedFrames = invspod(SPOD.P, A_reduced, SPOD.nDFT, SPOD.overlap);
    % Loop through each frame and add the mean field
    for k = 1:length(reconstructedFrames)
        reconstructedFrames(k, :, :) = reconstructedFrames(k, :, :) + reshape(SPOD.time_mean, [1, 2*h, w]);
    end



    % Create a VideoWriter object to write the video to a file
    videoFile = 'reconstructedFlow.avi';  % Video file name
    videoWriter = VideoWriter(videoFile, 'Uncompressed AVI');  % Create video writer object
    open(videoWriter);  % Open the video file for writing

    % Loop over all frames to create the video
    for k = 1:length(reconstructedFrames)
        
        % Get the current frame
        U = squeeze(reconstructedFrames(k, 1:h, :)); 
        V= squeeze(reconstructedFrames(k, h+1:end, :));
  
        % Display the current frame
        imagesc(U);
        axis equal tight;
        colorbar;

        % Write the current frame to the video
        writeVideo(videoWriter, getframe(gcf));  % Capture the frame and add to video

    end

    % Close the video file
    close(videoWriter);

    % Notify user
    disp(['Video saved as ', videoFile]);
end