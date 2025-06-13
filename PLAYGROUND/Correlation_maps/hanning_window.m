%% 1. Create a 16x16 Hanning window using the outer product of 1D Hanning windows
N = 16;
% Create a 1D Hanning window (MATLAB's hann function creates a symmetric Hanning window)
hanning1D = hann(N);

% Create the 2D Hanning window as the outer product
hanning2D = hanning1D * hanning1D';

% Plot the original 2D Hanning window
figure;
subplot(1,2,1);
imagesc(hanning2D);
colorbar;
title('2D Hanning Window');
axis equal tight;
xlabel('X'); ylabel('Y');

%% 2. Increase the Tapering Rate
% Raise the 2D window to a power greater than 1 to sharpen the drop-off
power = 2;  % Change this value (>1 gives a steeper taper)
hanning2D_sharp = hanning2D .^ power;

% Plot the sharpened window
subplot(1,2,2);
imagesc(hanning2D_sharp);
colorbar;
title(['2D Hanning Window, Raised to Power ' num2str(power)]);
axis equal tight;
xlabel('X'); ylabel('Y');

%% 3. Custom Radial Cosine Taper Centered at (9,9)
% Define coordinate grid
[x, y] = meshgrid(1:N, 1:N);
center = [9, 9];  % Center in MATLAB indexing

% Compute the Euclidean distance from the center
r = sqrt((x - center(1)).^2 + (y - center(2)).^2);

% Set the maximum radius at which the taper goes to zero.
% Here we choose the maximum distance from (9,9) to the array edge.
r_max = max([max(center-1), max(N - center)]);

% Create a radial cosine taper that decays to 0 at r_max
% The exponent p_radial controls the steepness of the taper.
p_radial = 1;  % Increase for a faster taper
taper = cos(pi * r / (2 * r_max)).^p_radial;
% Optionally force taper values to 0 beyond r_max if needed:
taper(r > r_max) = 0;

% Plot the radial cosine taper
figure;
imagesc(taper);
colorbar;
title('Custom Radial Cosine Taper Centered at (9,9)');
axis equal tight;
xlabel('X'); ylabel('Y');
