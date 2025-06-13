function b = build_b(node, dlet_c, dlet_n, dP_dx, dP_dy, f, h)
int_c = node.int_c;
west = node.west;   north = node.north;  south = node.south;  east = node.east;
ne = node.ne;       se = node.se;        nw = node.nw;        sw = node.sw;

% + load vector
b(int_c) = f(int_c)*h^2;
b(west(1,:))  = f(west(1,:))*h^2   + 2*h*dP_dx(west(1,:));
b(north(1,:)) = f(north(1,:))*h^2  - 2*h*dP_dy(north(1,:));
b(south(1,:)) = f(south(1,:))*h^2  + 2*h*dP_dy(south(1,:));
b(east(1,:))  = f(east(1,:))*h^2   - 2*h*dP_dx(east(1,:));
b(ne(1,:))    = 0.5*h^2*f(ne(1,:)) - h*dP_dx(ne(1,:)) - h*dP_dy(ne(1,:));
b(se(1,:))    = 0.5*h^2*f(se(1,:)) - h*dP_dx(se(1,:)) + h*dP_dy(se(1,:));
b(nw(1,:))    = 0.5*h^2*f(nw(1,:)) + h*dP_dx(nw(1,:)) - h*dP_dy(nw(1,:));
b(sw(1,:))    = 0.5*h^2*f(sw(1,:)) + h*dP_dx(sw(1,:)) + h*dP_dy(sw(1,:));

if ~isempty(dlet_c)
    b(dlet_c)     = -dlet_n(4,:);
end