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

samplingRate = [100,250,100,250,100,250,100,250,100,250]; % Hz 
dt = {3000*10^(-6), 1000*10^(-6), 3000*10^(-6), 1000*10^(-6), ...
      3000*10^(-6), 1000*10^(-6), 3000*10^(-6), 1000*10^(-6)};  % dt for different runs
min_expected_freq = 0.3;
imageCount = 18000;
caseImages = 3600;
endpoint = '';
x_loc = -5;
scaleFactor = 9.53;
CameraNo = 1;
close all;
i = 5;  % index to select a particular field from piv_result (for dimensions)
plot_type = 'v';
modes = [1,2];

for run = 1:length(base)
    % Load one calibration file to obtain image dimensions and mask
    VelData = load(fullfile(Calibrated, sprintf('%05d.mat', 1)));
    [h, w] = size(VelData.piv_result(i).ux);
    mask = VelData.piv_result(i).b_mask;
    row_average = mean(mask, 2);
    index = find(row_average > 0.8, 1);  % separation row index

    % Define locations for statistics and POD data
    statsloc = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                 ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated');
    dataloc_above = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                      ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated', 'Above');
    dataloc_below = fullfile(base{run}, 'Statistics', num2str(imageCount), ...
                      ['Cam', num2str(CameraNo)], 'Instantaneous', 'Calibrated', 'Below');

    % Load the POD statistics for the above and below regions
    Data_above = load(fullfile(dataloc_above, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));
    Data_below = load(fullfile(dataloc_below, ['POD_stats_', num2str(w), 'x', num2str(h), '.mat']));

    % (Optionally) adjust masks if needed
    mask_above = mask;
    mask_above(index-1:end, :) = 1;
    mask_below = mask;
    mask_below(1:index-1, :) = 1;
     
    % Load mean statistics (assumed to include mean_U and mean_V)
    StatsData = load(fullfile(statsloc, 'MeanStats16x16.mat'));
    mean_U = StatsData.mean_U;   % Expected to be vectorized for the region
    mean_V = StatsData.mean_V;
    u_prime_squared_calibrated = StatsData.U_prime_Uprime_mean;
    v_prime_squared_calibrated = StatsData.V_prime_Vprime_mean;
    TKE = 0.5 * (u_prime_squared_calibrated + v_prime_squared_calibrated);

    % Split TKE into the above and below regions
    TKE_masked_above = TKE(1:index-1, :);
    TKE_masked_below = TKE(index:end, :);
    Total_TKE_above = sum(TKE_masked_above(:));
    Total_TKE_below = sum(TKE_masked_below(:));

    % Compute uncertainty using the run-specific dt
    uncertainty = 0.1 / (scaleFactor/10^(-3) * dt{run});
    b_mask_cavity = mask(index-1:end, :);
    Area_below = sum(~b_mask_cavity(:));
    Area_above = sum(~mask_above(:));

    Noise_TKE_below = Area_below * uncertainty^2;
    Noise_TKE_above = Area_above * uncertainty^2;

    % Calculate signal-to-noise ratios (SNR)
    SNR_below = Total_TKE_below / Noise_TKE_below;
    SNR_above = Total_TKE_above / Noise_TKE_above;

    % Determine energy percentages and modes required
    Noise_percentage_below = 100 * (1 / SNR_below);
    Noise_percentage_above = 100 * (1 / SNR_above);
    real_signal_below = 100 - Noise_percentage_below;
    real_signal_above = 100 - Noise_percentage_above;
    k_below = find(Data_below.cumulativeEnergy > real_signal_below, 1, 'first');
    k_above = find(Data_above.cumulativeEnergy > real_signal_above, 1, 'first');
    
    disp(['The number of modes required to capture ', num2str(real_signal_below), '% of the energy is ', num2str(k_below), ' in cavity']);
    disp(['The number of modes required to capture ', num2str(real_signal_above), '% of the energy is ', num2str(k_above), ' above cavity']);
    
    % Preallocate arrays for the full field reconstructions
    
    for imNo = 1:imageCount
        %% Reconstruct the "above" (top) region
        % Reconstruct using the POD expansion
        UV_rec_above = Data_above.U_svd(:, 1:k_above) * ...
            (Data_above.S_svd(1:k_above, 1:k_above) * Data_above.V_svd(imNo, 1:k_above)');
        % Determine number of pixels in the above region
        n_above = (index - 1) * w;
        % Add the mean field; assume Data_above.mean_u and .mean_v are vectors of length n_above
        U_rec_above = UV_rec_above(1:n_above) + Data_above.mean_u;
        V_rec_above = UV_rec_above(n_above+1:end) + Data_above.mean_v;
        u_rec_above = reshape(U_rec_above, [index - 1, w]);
        v_rec_above = reshape(V_rec_above, [index - 1, w]);

        %% Reconstruct the "below" (cavity) region
        UV_rec_below = Data_below.U_svd(:, 1:k_below) * ...
            (Data_below.S_svd(1:k_below, 1:k_below) * Data_below.V_svd(imNo, 1:k_below)');
        % Determine number of pixels in the below region
        n_below = (h - index + 1) * w;
        U_rec_below = UV_rec_below(1:n_below) + Data_below.mean_u;
        V_rec_below = UV_rec_below(n_below+1:end) + Data_below.mean_v;
        u_rec_below = reshape(U_rec_below, [h - index + 1, w]);
        v_rec_below = reshape(V_rec_below, [h - index + 1, w]);

        %% Combine the two regions to form the full-field reconstruction
        u = [u_rec_above; u_rec_below];
        v = [v_rec_above; v_rec_below];


    end

    % (Optional) Save the reconstructed fields for this run
    % save(fullfile(base{run}, 'ReconstructedFlow.mat'), 'full_u', 'full_v', '-v7.3');
end
