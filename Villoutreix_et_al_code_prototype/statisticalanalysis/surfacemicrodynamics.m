function [m,s] = surfacemicrodynamics(sp,g,i,nts,L,selection, Options)
%DYNAMIQUEMICROVOLUME computes the averaged dynamics of the surface area 
%between consecutive mitoses for a given group of cell (population, 
%generation)

ngen6 = Options.ngen6;
% extract cells corresponding to the required population and generation
% with observed initial and final mitosis
indtemp = find(L{i}(:,2) == sp  & L{i}(:,5) == 1 & L{i}(:,6) == 1 & L{i}(:,7) == g);

m = [];
s = [];

% if the number of cells is at least greater than 95% of the expected 
% number of cells
if (length(indtemp)>=(0.95*2^(g-6)*ngen6(sp))),
    % storing the index of the cells
    nums = L{i}(indtemp,1);
    
    % initialization of the array where individual micro dynamics will be
    % stored
    surftemp1 = zeros(nts,length(nums));
    
    % initialization of the array where average and standard deviation over
    % the group of cells will be stored
    m = zeros(nts,1);
    s = zeros(nts,1);
    
    % for each cell
    for k = 1:length(nums),
        % store the complete micro dynamic of the cell surface area divided
        % by the averaged cell volume in the array voltemp
        surftemp = selection{i}(find(selection{i}(:,5) == nums(k)),7)/(L{i}(indtemp(k),10));
        % interpolate the dynamic of the volume on a common grid with nts
        % (default = 100) time steps
        surftemp1(:,k) = interp1((1:length(surftemp)),surftemp,(1:((length(surftemp)-1)/(nts-1)):length(surftemp)));
    end
    
    % for each time step of the common grid, compute the mean and standard
    % deviation of the normalized volume
    for t = 1:nts,
        m(t) = mean(surftemp1(t,find(surftemp1(t,:))));
        s(t) = std(surftemp1(t,find(surftemp1(t,:))));
    end
end
end