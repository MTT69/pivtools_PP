function p = intPolarPGrad(dpdx,dpdy,dx,theta)
% INTPOLARPGRAD - Obtains the integration path and integrated pressure
% field for the planar fields dpdx and dpdy for polar angle theta.
%
% Inputs:       dpdx - horizontal pressure gradient of size [ny nx]
%               dpdy - vertical pressure gradient of size [ny nx]
%                 dx - (uniform) grid spacing [meters]
%              theta - polar angle wrt horizontal (degrees)
%
% Output:          p - integrated pressure field

    % Begin by orienting the data
    if theta == 0 || theta == 360
        p = -fliplr(cumsum(fliplr(dpdx.*dx),2));
        return
    elseif theta == 90
        p = -cumsum(dpdy.*dx,1);
        return
    elseif theta == 180
        p = cumsum(dpdx.*dx,2); 
        return
    elseif theta == 270
        p = flipud(cumsum(flipud(dpdy.*dx),1));
        return
    % For non 90-degree angles, flip matrices to perform cumulative sum wrt
    % upper-left corner of the matrix
    elseif theta > 0 && theta < 90
        dpdx = fliplr(dpdx);
        dpdy = fliplr(dpdy);
        m = tand(theta);
    elseif theta > 90 && theta < 180
        dpdx = dpdx;
        dpdy = dpdy;
        m = tand(theta-90);
    elseif theta > 180 && theta < 270
        dpdx = flipud(dpdx);
        dpdy = flipud(dpdy);
        m = tand(theta-180);
    elseif theta > 270 && theta < 360
        dpdx = fliplr(flipud(dpdx));
        dpdy = fliplr(flipud(dpdy));
        m = tand(theta-270);
    end
    
    % Generate diagonal paths interpolated to measurement grid
    [ny,nx] = size(dpdx);
    nmx =  max([ny nx]);
    
    if m == 1
        
        diag_x = (0:nmx-1)+1;
        diag_y = (0:nmx-1)+1;
        uniq_y = diag_y;
        y_diff_v = ones(size(diag_x));
        x_diff_v = ones(size(diag_x));
        y_diff_h = ones(size(diag_x));
        x_diff_h = ones(size(diag_x));
        
    elseif m < 1 % x is independent variable
        
        i_x = 0:nmx-1;
        y = m.*(i_x);
        diag_y = round(y)+1;
        diag_x = i_x+1;
        [~,uniq_y] = unique(diag_y);
        uniq_y = [uniq_y; ones(length(diag_x)-length(uniq_y),1).*uniq_y(end)]';
        
        x_diff_v = diff([0 diag_x]);
        y_diff_v = diff([1 diag_y]);
        
        y_diff_h = diff([0 diag_x]);
        x_diff_h = diff([1 diag_x(uniq_y)]);
        
    elseif m > 1 % y is independent variable
        
        y = 0:nmx-1;
        i_x = (1/m).*(y);
        diag_x = round(i_x)+1;
        diag_y = y+1;
        [~,uniq_x] = unique(diag_x);
        uniq_x = [uniq_x; ones(length(diag_y)-length(uniq_x),1).*uniq_x(end)]';
        
        x_diff_v = diff([0 diag_y]);
        y_diff_v = diff([1 diag_y(uniq_x)]);
        
        y_diff_h = diff([0 diag_y]);
        x_diff_h = diff([1 diag_x]);
        
    end
    
    % Shift matrices according to diagonal summation path (vertical)
    % Account for non-singleton shifts for integration
    if m <= 1
        for i = 1:nx
            dpdx1_temp(:,i) = circshift(dpdx(:,i),-(diag_y(i)-1)).*x_diff_v(i);
            dpdy1_temp(:,i) = circshift(dpdy(:,i),-(diag_y(i)-1)).*y_diff_v(i);
        end
    elseif m > 1
        for i = 1:nx
            dpdx1_temp(:,i) = circshift(dpdx(:,i),-(diag_y(uniq_x(i))-1)).*x_diff_v(i);
            dpdy1_temp(:,i) = circshift(dpdy(:,i),-(diag_y(uniq_x(i))-1)).*y_diff_v(i);
        end
    end
    
    % Integrate (1 - lower domain triangle)
    % Correct sign depending on integration direction
    if theta > 0 && theta < 90
        px1 = -cumsum(dpdx1_temp.*dx,2);
        py1 = -cumsum(dpdy1_temp.*dx,2);
    elseif theta > 90 && theta < 180
        px1 = cumsum(dpdx1_temp.*dx,2);
        py1 = -cumsum(dpdy1_temp.*dx,2);
    elseif theta > 180 && theta < 270
        px1 = cumsum(dpdx1_temp.*dx,2);
        py1 = cumsum(dpdy1_temp.*dx,2);
    elseif theta > 270 && theta < 360
        px1 = -cumsum(dpdx1_temp.*dx,2);
        py1 = cumsum(dpdy1_temp.*dx,2);
    end
    
    % Unshift matrices
    if m <= 1
        for i = 1:nx
            px1(:,i) = circshift(px1(:,i),(diag_y(i)-1));
            py1(:,i) = circshift(py1(:,i),(diag_y(i)-1));
        end
    elseif m > 1
        for i = 1:nx
            px1(:,i) = circshift(px1(:,i),(diag_y(uniq_x(i))-1));
            py1(:,i) = circshift(py1(:,i),(diag_y(uniq_x(i))-1));
        end
    end
    
    % Shift matrices according to diagonal summation path (horizontal)
    % Account for non-singleton shifts for integration
    for i = 1:ny
        if m <= 1
            dpdx2_temp(i,:) = circshift(dpdx(i,:),[0 -(diag_x(uniq_y(i))-1)]).*x_diff_h(i);
            dpdy2_temp(i,:) = circshift(dpdy(i,:),[0 -(diag_x(uniq_y(i))-1)]).*y_diff_h(i);
        elseif m > 1
            dpdx2_temp(i,:) = circshift(dpdx(i,:),[0 -(diag_x(i)-1)]).*x_diff_h(i);
            dpdy2_temp(i,:) = circshift(dpdy(i,:),[0 -(diag_x(i)-1)]).*y_diff_h(i);
        end
    end
    
    % Integrate (2 - upper domain triangle)
    % Correct sign depending on integration direction
    if theta > 0 && theta < 90
        px2 = -cumsum(dpdx2_temp.*dx,1);
        py2 = -cumsum(dpdy2_temp.*dx,1);
    elseif theta > 90 && theta < 180
        px2 = cumsum(dpdx2_temp.*dx,1);
        py2 = -cumsum(dpdy2_temp.*dx,1);
    elseif theta > 180 && theta < 270
        px2 = cumsum(dpdx2_temp.*dx,1);
        py2 = cumsum(dpdy2_temp.*dx,1);
    elseif theta > 270 && theta < 360
        px2 = -cumsum(dpdx2_temp.*dx,1);
        py2 = cumsum(dpdy2_temp.*dx,1);        
    end   
    
    % Unshift matrices
    for i = 1:ny
        if m <= 1
            px2(i,:) = circshift(px2(i,:),[0 (diag_x(uniq_y(i))-1)]);
            py2(i,:) = circshift(py2(i,:),[0 (diag_x(uniq_y(i))-1)]);
        elseif m > 1
            px2(i,:) = circshift(px2(i,:),[0 (diag_x(i)-1)]);
            py2(i,:) = circshift(py2(i,:),[0 (diag_x(i)-1)]);
        end
    end
    
    % Extract triangles
    if m == 1
        px1 = tril(px1,-1);
        px2 = triu(px2);
        py1 = tril(py1,-1);
        py2 = triu(py2);
    elseif m > 1 
        lower = zeros(ny,nx);
        for i = 1:ny
            lower(i,1:diag_x(i)-1) = 1;
        end
        lower = lower(1:ny,1:nx);
        upper = imcomplement(lower);
        px1 = px1.*lower;
        px2 = px2.*upper;
        py1 = py1.*lower;
        py2 = py2.*upper;
    elseif m < 1 
        upper = zeros(ny,nx);
        for i = 1:nx
           upper(1:diag_y(i)-1,i) = 1; 
        end
        upper = upper(1:ny,1:nx);
        lower = imcomplement(upper);
        px1 = px1.*lower;
        px2 = px2.*upper;
        py1 = py1.*lower;
        py2 = py2.*upper;
    end
    
    % Combine, reorient, and stitch paths
    if theta > 0 && theta < 90
        px = fliplr(sum(cat(3,px1,px2),3));
        py = fliplr(sum(cat(3,py1,py2),3));
        p = px+circshift(py,[0 1]);
    elseif theta > 90 && theta < 180
        px = sum(cat(3,px1,px2),3);
        py = sum(cat(3,py1,py2),3);
        p = px+circshift(py,[0 -1]);
    elseif theta > 180 && theta < 270
        px = flipud(sum(cat(3,px1,px2),3));
        py = flipud(sum(cat(3,py1,py2),3));
        p = px+circshift(py,[0 -1]);
    elseif theta > 270 && theta < 360
        px = fliplr(flipud(sum(cat(3,px1,px2),3)));
        py = fliplr(flipud(sum(cat(3,py1,py2),3)));
        p = px+circshift(py,[0 1]);
    end
    
end
