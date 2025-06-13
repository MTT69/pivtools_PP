function [dP_dx, dP_dy, f] = TH(X,Y,U_MEAN,V_MEAN,U,V,rho,nu,node)
up = U - U_MEAN;
vp = V - V_MEAN;

% first order derivatives
[U_dx, U_dy] = grad2D(U,X,Y,node);    
[V_dx, V_dy] = grad2D(V,X,Y,node);   
[up_dx, up_dy] = grad2D(up,X,Y);
[vp_dx, vp_dy] = grad2D(vp,X,Y);

% second order derivatives
[U_dx2, ~] =  grad2D(U_dx,X,Y,node);     
[~, U_dy2] = grad2D(U_dy,X,Y,node);      
[V_dx2, ~] =  grad2D(V_dx,X,Y,node);     
[~, V_dy2] = grad2D(V_dy,X,Y,node);    

% pressure gradient TH1
dP_dx = -rho*(-(U_MEAN.*up_dx+V_MEAN.*up_dy) + (U.*U_dx+V.*U_dy) - nu*(U_dx2+U_dy2));
dP_dy = -rho*(-(U_MEAN.*vp_dx+V_MEAN.*vp_dy) + (U.*V_dx+V.*V_dy) - nu*(V_dx2+V_dy2));


% poisson equation for pressure
[gx, ~] = grad2D(dP_dx,X,Y,node);
[~ , gy] = grad2D(dP_dy,X,Y,node);
f = (gx + gy);