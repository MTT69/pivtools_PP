%% assumes figure has been opened
% h=gcf;
% ax = findobj(h, 'Type', 'Axes');  % Find all axes in the figure
% lower_limit = 0;
% upper_limit = 0.02;

% title(ax(2), [''], 'FontSize', 20, 'Interpreter', 'latex');


% set(ax(2), 'FontSize', 18);
% set(ax(2), 'XTick', [], 'YTick', []);
% colorbar(ax(2),'FontSize', 18);
% cb = colorbar(ax(2), 'FontSize', 18);

% 
% set(ax(1), 'CLim', [lower_limit upper_limit]);
% custommap = redbluezero(lower_limit,upper_limit);
% 
% colormap(ax(1), custommap);  % You can choose any other colormap like 'parula', 'hot', 'cool', etc.

% 
axesHandles = findall(gcf, 'Type', 'axes');
for i = 1:length(axesHandles)
    xlim(axesHandles(i), [0, 5]);

end






% Assuming ax is your axis handle
% xLimits = xlim(ax(2));  % Get the current x-axis limits
% 
% % Round the limits to the nearest integers
% minX = ceil(xLimits(1));  % Round up the lower limit
% maxX = floor(xLimits(2));  % Round down the upper limit

% % Set the x-ticks to increment by 1 from the rounded limits
% % xticks(ax(2), minX:1:maxX);
% yticks = y_co_ords_large(:,1);
% set(ax(2), 'YTick', yticks(1:10:end));
% daspect([1 1 1])