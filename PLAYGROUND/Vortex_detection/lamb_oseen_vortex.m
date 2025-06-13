%% Synthetic Vortex Flow and Gamma2 Vortex Centre Identification
close all; clear; clc;

%% Parameters for the synthetic vortex flow
Gamma = 1.0;    % Circulation [m^2/s]
nu = 0.01;      % Kinematic viscosity [m^2/s]
t = 1.0;        % Time [s]

% Meshgrid
x_min = -1; x_max = 1;  % X-axis limits
y_min = -1; y_max = 1;  % Y-axis limits
n_points = 100;         % Number of points in each dimension
[x, y] = meshgrid(linspace(x_min, x_max, n_points), linspace(y_min, y_max, n_points));
y = flipud(y);  % flip y to match typical plotting orientation

% Radial distance from the vortex centre (assumed at 0,0)
r = sqrt(x.^2 + y.^2);
r(r==0) = 1e-10;  % avoid division by zero

% Velocity components (Lamb–Oseen vortex)
v_theta = (Gamma ./ (2 * pi * r)) .* (1 - exp(-r.^2 / (4 * nu * t)));
u = -v_theta .* (y ./ r);  % x-component
v =  v_theta .* (x ./ r);   % y-component

%% (Optional) Compute vorticity for reference (not used for vortex centre detection)
[dx,~] = gradient(x);
[~,dy] = gradient(y);
dy = -dy;  % because y was flipped

[dvx, ~] = gradient(v);
[~, duy] = gradient(u);
duy = -duy; % due to y orientation

vorticity = (dvx ./ dx - duy ./ dy);

%% Compute the Gamma2 Field
% The Gamma2 criterion (Graftieaux et al., 2001) identifies vortex centres as locations
% where the local swirling strength (computed over a neighborhood) is near +1 (or -1).
%
% For each candidate point (x0,y0), we compute:
%
%   Gamma2(x0,y0) = (1/N) * sum_{neighbors in S} { ( (r x v) / (|r| |v|) ) }
%
% where r = (x - x0, y - y0) is the relative position vector,
%       v = (u, v) is the local velocity,
% and the sum is over all points in a circular neighborhood S (of radius R).
%
% We now loop over the grid and compute the Gamma2 value at each grid point.
%
R = 0.1;  % Radius for the neighborhood (adjust as needed)
gamma2Field = zeros(size(x));  % Preallocate

% Loop over each grid point
for i = 1:n_points
    for j = 1:n_points
        x0 = x(i,j);
        y0 = y(i,j);
        gamma2Field(i,j) = computeGamma2(x, y, u, v, x0, y0, R);
    end
end

%% Identify Vortex Centres
% In an ideal vortex, Gamma2 is near +1 (or –1) at the vortex centre.
% Here we use a threshold (e.g., abs(Gamma2) > 0.9) and also look for local maxima.
threshold = 0.9;
% Use MATLAB's imregionalmax to find local maxima in the absolute Gamma2 field.
BW = imregionalmax(abs(gamma2Field));
% Apply the threshold to filter out spurious peaks.
BW = BW & (abs(gamma2Field) > threshold);

% Extract the positions of the detected vortex centres.
vortexCentersX = x(BW);
vortexCentersY = y(BW);

%% Plot Results
figure;
subplot(2,2,1);
contourf(x, y, u, 20); colorbar;
title('u-velocity'); xlabel('x'); ylabel('y');

subplot(2,2,2);
contourf(x, y, v, 20); colorbar;
title('v-velocity'); xlabel('x'); ylabel('y');

subplot(2,2,3);
contourf(x, y, vorticity, 20); colorbar;
title('Vorticity'); xlabel('x'); ylabel('y');

subplot(2,2,4);
contourf(x, y, gamma2Field, 20); colorbar; hold on;
plot(vortexCentersX, vortexCentersY, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
title('Gamma2 Field and Detected Vortex Centres');
xlabel('x'); ylabel('y'); hold off;



%% Helper Function: computeGamma2
function gamma2 = computeGamma2(x, y, u, v, x0, y0, R)
    % Create a mask for points within a circle of radius R centered at (x0, y0)
    mask = (x - x0).^2 + (y - y0).^2 <= R^2;
    
    % Relative position vectors (r = [rx, ry])
    rx = x(mask) - x0;
    ry = y(mask) - y0;
    
    % Velocity components in the neighborhood
    u_local = u(mask);
    v_local = v(mask);
    
    % Compute the norms of the relative positions and local velocities.
    % Add eps to avoid division by zero (especially at the centre point)
    r_norm = sqrt(rx.^2 + ry.^2) + eps;
    v_norm = sqrt(u_local.^2 + v_local.^2) + eps;
    
    % Compute the (scalar) cross product (in 2D, a x b = a_x*b_y - a_y*b_x)
    cross_prod = rx .* v_local - ry .* u_local;
    
    % The normalized contribution from each point in the neighborhood.
    gamma_vals = cross_prod ./ (r_norm .* v_norm);
    
    % The Gamma2 value is the mean of these contributions.
    gamma2 = mean(gamma_vals);
end
