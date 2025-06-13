function [dP_dx, dP_dy, f] = RANS(X,Y,U_MEAN,V_MEAN,uv_MEAN,U_VAR,V_VAR,rho,nu,node)
% first order derivatives
[U_MEAN_dx, U_MEAN_dy] = grad2D(U_MEAN,X,Y,node);  
[V_MEAN_dx, V_MEAN_dy] = grad2D(V_MEAN,X,Y,node);  
[uv_dx, uv_dy] = grad2D(uv_MEAN,X,Y,node);  
[UVAR_dx, ~] = grad2D(U_VAR,X,Y,node);   
[~, VVAR_dy] = grad2D(V_VAR,X,Y,node);       

% second order derivatives
[U_MEAN_dx2, ~] =  grad2D(U_MEAN_dx,X,Y,node);     
[~, U_MEAN_dy2] = grad2D(U_MEAN_dy,X,Y,node);      
[V_MEAN_dx2, ~] =  grad2D(V_MEAN_dx,X,Y,node);    
[~, V_MEAN_dy2] = grad2D(V_MEAN_dy,X,Y,node);
% 
% figure; 
% subplot(2,2,1); imagesc(U_MEAN_dx2);
% subplot(2,2,2); imagesc(U_MEAN_dy2);
% subplot(2,2,3); imagesc(V_MEAN_dx2);
% subplot(2,2,4); imagesc(V_MEAN_dy2);

% pressure gradient
dP_dx = -rho*(U_MEAN.*U_MEAN_dx + V_MEAN.*U_MEAN_dy + UVAR_dx + uv_dy - nu*(U_MEAN_dx2+U_MEAN_dy2));
dP_dy = -rho*(U_MEAN.*V_MEAN_dx + V_MEAN.*V_MEAN_dy + VVAR_dy + uv_dx - nu*(V_MEAN_dx2+V_MEAN_dy2) );

% poisson equation for pressure
[gx, ~] = grad2D(dP_dx,X,Y,node);
[~ , gy] = grad2D(dP_dy,X,Y,node);
f = (gx + gy);



