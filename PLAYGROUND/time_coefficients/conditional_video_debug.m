figure;
x =x(end-72:end,:);
y=y(end-72:end,:);
mask_below = mask_below(end-72:end,:);
for j = 1:243
    disp(j)
    % Process each matrix and apply the mask
    A = developmentVideo_bottom_mean{j};
    A(mask_below == 1) = 0;
    
    B = developmentVideo_ux_bottom_mean{j};
    B(mask_below == 1) = 0;
    
    C = developmentVideo_uy_bottom_mean{j};
    C(mask_below == 1) = 0;
    
    D = gamma1(x, y, B, C, 10);
    D(mask_below == 1) = 0;

    % Plot in a 2x2 layout
    clf; % Clear figure for next frame
    
    subplot(2, 2, 1);
    imagesc(A, [-0.03, 0.03]);
    colorbar;
    title('Development Video Bottom Mean');
    
    subplot(2, 2, 2);
    imagesc(B);
    colorbar;
    title('Development Video Ux Bottom');
    
    subplot(2, 2, 3);
    imagesc(C);
    colorbar;
    title('Development Video Uy Bottom');
    
    subplot(2, 2, 4);
    imagesc(D);
    colorbar;
    title('Gamma1 Function Output');
    pause(0.05)


    

end


for j = 1:243
    disp(j)
    % Process each matrix and apply the mask
    A = developmentVideo_top_mean{j};

    
    B = developmentVideo_ux_top_mean{j};

    
    C = developmentVideo_uy_top_mean{j};

   

    % Plot in a 2x2 layout
    clf; % Clear figure for next frame
    
    subplot(2, 2, 1);
    imagesc(A, [-0.03, 0.03]);
    colorbar;
    title('Development Video Bottom Mean');
    
    subplot(2, 2, 2);
    imagesc(B);
    colorbar;
    title('Development Video Ux Bottom');
    
    subplot(2, 2, 3);
    imagesc(C);
    colorbar;
    title('Development Video Uy Bottom');
 
    pause(0.05)

    

end




mask_below = mask_below(end-72:end,:);