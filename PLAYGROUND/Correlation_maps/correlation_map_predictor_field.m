% Fix typo in variable name
run = 7;
predx = piv_result(run).Predictor_Field(:,:,2);
predy = piv_result(run).Predictor_Field(:,:,1); % Corrected from Predcitor_Field
ux = piv_result(run).ux;
uy=piv_result(run).uy;

% Extract only the bottom 10 rows
bottom_rows = size(predx, 1) - 10:size(predx, 1);
bottom_cols = ceil(size(predx,2)/2)-25:ceil(size(predx,2)/2)+25;
predx = predx(bottom_rows, bottom_cols);
predy = predy(bottom_rows, bottom_cols);
ux = ux(bottom_rows, bottom_cols);
uy = uy(bottom_rows, bottom_cols);

% Plot for x-component
figure(1);

% Calculate common color limits for x-component
x_min = min([predx(:); ux(:)]);
x_max = max([predx(:); ux(:)]);

subplot(3,1,1);
imagesc(predx);
colorbar;
caxis([x_min, x_max]);
title('Predictor Field - X Component (Bottom 10 vectors)');
axis equal tight;

subplot(3,1,2);
imagesc(ux);
colorbar;
caxis([x_min, x_max]);
title('Velocity Field - X Component (Bottom 10 vectors)');
axis equal tight;

subplot(3,1,3);
imagesc(ux - predx);
colorbar;
caxis([-1, 1]);  % Set comparison range between -1 and 1
title('Difference (ux - predx) (Bottom 10 vectors)');
axis equal tight;

% Plot for y-component
figure(2);

% Calculate common color limits for y-component
y_min = min([predy(:); uy(:)]);
y_max = max([predy(:); uy(:)]);

subplot(3,1,1);
imagesc(predy);
colorbar;
caxis([y_min, y_max]);
title('Predictor Field - Y Component (Bottom 10 vectors)');
axis equal tight;

subplot(3,1,2);
imagesc(uy);
colorbar;
caxis([y_min, y_max]);
title('Velocity Field - Y Component (Bottom 10 vectors)');
axis equal tight;

subplot(3,1,3);
imagesc(uy - predy);
colorbar;
caxis([-0.5, 0.5]);  % Set comparison range between -0.5 and 0.5
title('Difference (uy - predy) (Bottom 10 vectors)');
axis equal tight;
