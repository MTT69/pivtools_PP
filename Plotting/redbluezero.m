function c = redbluezero(lowerlim, upperlim)
% This function generates a colormap that transitions from shades of blue to 
% white and then to shades of red, with the midpoint at zero (the color white). 
% The colormap adjusts the gradients to respect asymmetric input ranges, ensuring
% that the visual transition matches the actual data range while keeping white at zero.
%
% Inputs:
%   lowerlim   - The lower limit of the data range (negative value).
%   upperlim   - The upper limit of the data range (positive value).
%
% Outputs:
%   c          - An M-by-3 matrix representing the RGB values of the colormap.
%                Transitions are asymmetric to reflect the data range, but white remains at zero.

% Define the number of colors in the colormap
m = 256;

% Compute proportions for negative and positive ranges
total_range = abs(lowerlim) + upperlim;
negative_proportion = abs(lowerlim) / total_range;
positive_proportion = upperlim / total_range;

% Calculate the number of colors for each range
m_neg = floor(m * negative_proportion);
m_pos = m - m_neg;

% Blue to white for the negative range
r_neg = linspace(0, 1, m_neg)';
g_neg = linspace(0, 1, m_neg)';
b_neg = ones(m_neg, 1);

% White to red for the positive range
r_pos = ones(m_pos, 1);
g_pos = linspace(1, 0, m_pos)';
b_pos = linspace(1, 0, m_pos)';

% Combine the two parts to create the full colormap
r = [r_neg; r_pos];
g = [g_neg; g_pos];
b = [b_neg; b_pos];

% Combine RGB channels
c = [r g b];
end
