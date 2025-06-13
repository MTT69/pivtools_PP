function [M_dx, M_dy] = grad2D(M, X, Y, node)
if nargin == 3 % if nodes have not been identified yet
    pad = 1;
    M = padarray(M,[pad,pad],nan,'both');
    node = nodeID(M);
end

int_c = node.int_c; int_n = node.int_n;
west = node.west;   north = node.north;  south = node.south;  east = node.east;
ne = node.ne;       se = node.se;        nw = node.nw;        sw = node.sw;
M_dx = nan(size(M));
M_dy = nan(size(M));
dx = (X(2)-X(1));
dy = (Y(2)-Y(1));

% central differences on internal nodes
M_dx(int_c) = (M(int_n(4,:)) - M(int_n(1,:)))/(2*dx);
M_dy(int_c) = (M(int_n(3,:)) - M(int_n(2,:)))/(2*dy);

% forward/backward + central differences on nodes at plane surface
M_dx(west(1,:)) = (M(west(2,:)) - M(west(1,:)))/dx;
M_dy(west(1,:)) = (M(west(4,:)) - M(west(3,:)))/(2*dy);

M_dx(north(1,:)) = (M(north(4,:)) - M(north(3,:)))/(2*dx);
M_dy(north(1,:)) = (M(north(2,:)) - M(north(1,:)))/dy;

M_dx(south(1,:)) = (M(south(4,:)) - M(south(3,:)))/(2*dx);
M_dy(south(1,:)) = (M(south(1,:)) - M(south(2,:)))/dy;

M_dx(east(1,:)) = (M(east(1,:)) - M(east(2,:)))/dx;
M_dy(east(1,:)) = (M(east(4,:)) - M(east(3,:)))/(2*dy);

% forward and backward derivatives on external corner nodes
M_dx(nw(1,:)) = (M(nw(3,:)) - M(nw(1,:)))/dx; 
M_dy(nw(1,:)) = (M(nw(2,:)) - M(nw(1,:)))/dy;
M_dx(sw(1,:)) = (M(sw(3,:)) - M(sw(1,:)))/dx; 
M_dy(sw(1,:)) = (M(sw(1,:)) - M(sw(2,:)))/dy;
M_dx(ne(1,:)) = (M(ne(1,:)) - M(ne(2,:)))/dx; 
M_dy(ne(1,:)) = (M(ne(3,:)) - M(ne(1,:)))/dy;
M_dx(se(1,:)) = (M(se(1,:)) - M(se(2,:)))/dx; 
M_dy(se(1,:)) = (M(se(1,:)) - M(se(3,:)))/dy;

if nargin == 3
    M_dx = M_dx(pad+1:end-pad,pad+1:end-pad);
    M_dy = M_dy(pad+1:end-pad,pad+1:end-pad);
end
