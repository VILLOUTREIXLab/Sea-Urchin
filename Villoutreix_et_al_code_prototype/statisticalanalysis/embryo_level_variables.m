function [ Nbcells, Volcells, Surfcells ] = embryo_level_variables( Data, Options)
%EMBRYO_LEVEL_VARIABLES computes embryo level variables:
%total number of cells, total cell volume, total cell surface area

% Number of cells - cells disappearing from the cell tracking are maintained
% in the total count since there is no cell death in the considered period
Nbcells = zeros(max(Options.tmax),5,Options.Nb);

for i = 1:Options.Nb,%for each embryo
    for k = 1:length(Data.L{i}(:,1)),% going through all cells
        % total number of cells is incremented by thourghout the complete
        % cell cycle
        % for the whole embryo
        Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),5,i) = 1+Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),5,i);
        % for each population
        Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),Data.L{i}(k,2),i) = 1+Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),Data.L{i}(k,2),i);

    end
end

% volume and surface
Volcells = zeros(max(Options.tmax),5,Options.Nb);
Surfcells = zeros(max(Options.tmax),5,Options.Nb);

for i = 1:Options.Nb,% for each embryo
    for t = 1:Options.tmax(i), % for each time step
        % the total cell volume is computed as the sum of each individual
        % cell volume alive at time t
        % for the whole embryo
        indtemp = find(Data.selection{i}(:,3) == t);
        Volcells(t,5,i) = sum(Data.selection{i}(indtemp,6));
        Surfcells(t,5,i) = sum(Data.selection{i}(indtemp,7));
        % and for each population
        for j = 1:4,
            indtemp = find(Data.selection{i}(:,3) == t & Data.selection{i}(:,2) == j);
            Volcells(t,j,i) = sum(Data.selection{i}(indtemp,6));
            Surfcells(t,j,i) = sum(Data.selection{i}(indtemp,7));
        end
    end
end

end

