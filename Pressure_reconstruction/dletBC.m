function [dlet_c, dlet_n, node] = dletBC(dlet_ind, dlet_d, dlet_value, node)
    dlet_c = [];
    dlet_n = [];
for i = 1:length(dlet_ind)
    switch dlet_d{i}
        case 'west'
            dlet_c = [dlet_c,dlet_ind{i}];
            [~,aux] = intersect(node.west(1,:), dlet_ind{i});
            for j = 1:length(dlet_ind{i})
                dlet_n = [dlet_n, [node.west(2:4,aux(j)); dlet_value{i}(j)] ];
            end
            node.west(:,aux) = [];
        case 'north'
            dlet_c = [dlet_c,dlet_ind{i}];
            [~,aux] = intersect(node.north(1,:), dlet_ind{i});
            for j = 1:length(dlet_ind{i})
                dlet_n = [dlet_n, [node.north(2:4,aux(j)); dlet_value{i}(j)] ];
            end
            node.north(:,aux) = [];
        case 'south'
            dlet_c = [dlet_c,dlet_ind{i}];
            [~,aux] = intersect(node.south(1,:), dlet_ind{i});
            for j = 1:length(dlet_ind{i})
                dlet_n = [dlet_n, [node.south(2:4,aux(j)); dlet_value{i}(j)] ];
            end
            node.south(:,aux) = [];
        case 'east'
            dlet_c = [dlet_c,dlet_ind{i}];
            [~,aux] = intersect(node.east(1,:), dlet_ind{i});
            for j = 1:length(dlet_ind{i})
                dlet_n = [dlet_n, [node.east(2:4,aux(j)); dlet_value{i}(j)] ];
            end
            node.east(:,aux) = [];
    end
end