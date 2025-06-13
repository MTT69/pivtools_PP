function [I,J,V] = build_A(node, dlet_c, dlet_n)
int_c = node.int_c; int_n = node.int_n;
west = node.west;   north = node.north;  south = node.south;  east = node.east;
ne = node.ne;       se = node.se;        nw = node.nw;        sw = node.sw;
    
% + sparse matrix
aux1  = padarray(int_c,3,'symmetric','pre');
aux1  = reshape(aux1,[numel(aux1),1]);
int_n = reshape(int_n,[numel(aux1),1]);

%         int nodes,        int nbours,
sparse = {[int_c', int_c', -4*ones(size(int_c,2),1)],       [aux1, int_n, ones(size(int_n,1),1)],                   [],                                                 [];...
    % west nodes,               east nbour,                 remaining nbours,
    [west(1,:)',west(1,:)', -4*ones(size(west,2),1)],       [west(1,:)',west(2,:)',     2*ones(size(west,2),1)],    [west(1,:)',west(3,:)', ones(size(west,2),1)],      [west(1,:)',west(4,:)', ones(size(west,2),1)];...
    % north nodes,              south nbour,                remaining nbours,
    [north(1,:)',north(1,:)', -4*ones(size(north,2),1)],	[north(1,:)',north(2,:)',   2*ones(size(north,2),1)],	[north(1,:)',north(3,:)', ones(size(north,2),1)],	[north(1,:)',north(4,:)', ones(size(north,2),1)];...
    % south nodes,              north nbour,                remaining nbours,
    [south(1,:)',south(1,:)', -4*ones(size(south,2),1)],	[south(1,:)',south(2,:)',   2*ones(size(south,2),1)],	[south(1,:)',south(3,:)', ones(size(south,2),1)],	[south(1,:)',south(4,:)', ones(size(south,2),1)];...
    % east nodes,               west nbour,                 remaining nbours,
    [east(1,:)',east(1,:)', -4*ones(size(east,2),1)],       [east(1,:)',east(2,:)',     2*ones(size(east,2),1)],	[east(1,:)',east(3,:)', ones(size(east,2),1)],      [east(1,:)',east(4,:)', ones(size(east,2),1)]...
    % ne nodes,             nbours
    [ne(1,:)', ne(1,:)', -2*ones(size(ne,2),1)],            [ne(1,:)', ne(2,:)', ones(size(ne,2),1)],               [ne(1,:)', ne(3,:)', ones(size(ne,2),1)],           [];...
    % se nodes,             nbours
    [se(1,:)', se(1,:)', -2*ones(size(se,2),1)],            [se(1,:)', se(2,:)', ones(size(se,2),1)],               [se(1,:)', se(3,:)', ones(size(se,2),1)],           [];...
    % nw nodes,             nbours
    [nw(1,:)', nw(1,:)', -2*ones(size(nw,2),1)],            [nw(1,:)', nw(2,:)', ones(size(nw,2),1)],               [nw(1,:)', nw(3,:)', ones(size(nw,2),1)],           [];...
    % sw nodes,             nbours
    [sw(1,:)', sw(1,:)', -2*ones(size(sw,2),1)],            [sw(1,:)', sw(2,:)', ones(size(sw,2),1)],               [sw(1,:)', sw(3,:)', ones(size(sw,2),1)],           []};

if ~isempty(dlet_c)
    aux2  = padarray(dlet_c,2,'symmetric','pre');
    aux2  = reshape(aux2,[numel(aux2),1]);
    dlet_n = reshape(dlet_n(1:3,:),[numel(aux2),1]);
    % nodes with Dirichlet BC
    sparse(end+1,:) =  {[dlet_c', dlet_c', -4*ones(size(dlet_c,2),1)], [aux2, dlet_n, ones(size(dlet_n,1),1)], [], []};
end

I = []; J = []; V = [];
for i = 1:numel(sparse)
    if isempty(sparse{i}) == 0
        I = [I;sparse{i}(:,1)];
        J = [J;sparse{i}(:,2)];
        V = [V;sparse{i}(:,3)];
    end
end