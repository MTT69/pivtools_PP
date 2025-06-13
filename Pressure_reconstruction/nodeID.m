function node = nodeID(dummy)

% + find neighbouring nodes
[r,c] = size(dummy);
offset = [-r, -1, +1, +r];                   % vector offset
idx_c  = 1:r*c;                              % center nodes
idx_n  = bsxfun(@plus, (1:r*c), offset(:));  % linear indices of
% neighbouring nodes [east, north, south, west]

% + remove center nan nodes and their respective neighbours
idx_c(ind2sub(size(dummy),find(isnan(dummy)==1))) = [];
idx_n(:,ind2sub(size(dummy),find(isnan(dummy)==1))) = [];

% + remove neighbouring nan nodes
idx_n(isnan(dummy(idx_n(:))) == 1) = nan;

% + identify type of nodes and its neighbours
count  = nansum(idx_n./idx_n,1);
node.int_c  = idx_c(:,count == 4);     % interior nodes
node.int_n  = idx_n(:,count == 4);
bnd_c1 = idx_c(:,count == 3);     % nodes at plane surfaces
bnd_n1 = idx_n(:,count == 3);
bnd_c2 = idx_c(:,count == 2);     % nodes at external corner
bnd_n2 = idx_n(:,count == 2);

% + breakdown nodes at plane surface
node.west  = [bnd_c1(isnan(bnd_n1(1,:)));	bnd_n1(4,isnan(bnd_n1(1,:)));	bnd_n1([2,3],isnan(bnd_n1(1,:)))];
node.north = [bnd_c1(isnan(bnd_n1(2,:)));	bnd_n1(3,isnan(bnd_n1(2,:)));	bnd_n1([1,4],isnan(bnd_n1(2,:)))];
node.south = [bnd_c1(isnan(bnd_n1(3,:)));	bnd_n1(2,isnan(bnd_n1(3,:)));	bnd_n1([1,4],isnan(bnd_n1(3,:)))];
node.east  = [bnd_c1(isnan(bnd_n1(4,:)));	bnd_n1(1,isnan(bnd_n1(4,:)));	bnd_n1([2,3],isnan(bnd_n1(4,:)))];

% + breakdown nodes at external corners
aux2 = isnan(bnd_n2(1,:)) & isnan(bnd_n2(2,:));
node.nw = [bnd_c2(aux2);	bnd_n2([3,4],aux2)];
aux2 = isnan(bnd_n2(1,:)) & isnan(bnd_n2(3,:));
node.sw = [bnd_c2(aux2);	bnd_n2([2,4],aux2)];
aux2 = isnan(bnd_n2(2,:)) & isnan(bnd_n2(4,:));
node.ne = [bnd_c2(aux2);	bnd_n2([1,3],aux2)];
aux2 = isnan(bnd_n2(3,:)) & isnan(bnd_n2(4,:));
node.se = [bnd_c2(aux2);	bnd_n2([1,2],aux2)];

% field  = im2bw(abs(dummy), 0);
% figure(1);
% imshow(field); hold on
% [a,b] = ind2sub([r,c],node.int_c);
% scatter(b,a,7,'r','filled')
% [a,b] = ind2sub([r,c],bnd_c1);
% scatter(b,a,7,'b','filled')
% [a,b] = ind2sub([r,c],bnd_c2);
% scatter(b,a,7,'g','filled')
% axis equal

