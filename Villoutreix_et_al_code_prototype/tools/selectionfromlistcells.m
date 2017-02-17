function [selection] = selectionfromlistcells(list_cells, ind_volume, ind_surface, time)
%SELECTIONFROMLISTCELLS generates a selection file from an L-array (list
%of global cell features)

%the columns of selection correspond for each
%cell at each time step to the following features: 1) cell population 2)
% generation 3) id mother, 5) volume, 6) surface area
selection = zeros(1,7);

for k = 1:length(time),
    indtemp = find(list_cells(:,5)<=time(k) & list_cells(:,6)>time(k));
    selection = [selection;[indtemp,k*ones(length(indtemp),1),list_cells(indtemp,[1,2,3,ind_volume,ind_surface])]];
end

end