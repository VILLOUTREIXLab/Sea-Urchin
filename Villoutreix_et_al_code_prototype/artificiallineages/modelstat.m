function [ nbcellsstatmean,volcellsstatmean,surfcellsstatmean ] = modelstat( i,ng,ndraw,Data, Options, pop_parameters)
%MODELSTAT computes artificial cell lineages based on population level
%statistics and cell volume and surface area microdynamics - it returns
%embryo-level dynamics such as the total number of cells, the total cell
%volume and the total cell surface

% initialization of useful constants and variables
tinit = Data.rescaledtime{i}(1)/60;
tfin = Data.rescaledtime{i}(Options.tmax(i))/60;
% time grid
time = tinit:tfin;

g0 = Options.g0(i,:);

% initial number of cells for each population
n = Options.ngen6.*2.^(g0 - 6*ones(1,4));

% initialization of the arrays that will keep track of the embryo-level
% dynamics from the random cell lineages
nbcellsstat = zeros(length(time),5,ndraw);
volcellsstat = zeros(length(time),5,ndraw);
surfcellsstat = zeros(length(time),5,ndraw);



% we begin by computing the volume and surface area micro dynamics that
% will be used later when generating artificial cell lineages

%-------------
% Computing volume and surface area micro dynamics
%-------------
% number of time step of the grid on which the micro dynamics will be
% computed
nts = 100;
vol_micro_dyn = zeros(nts,4,10);
surf_micro_dyn = zeros(nts,4,10);
% loop on the populations
for j = 1:4,
    % loop on the generations
    for g = 6:10,
        % Computing averaged micro dynamics volume
        [m,~] = volumemicrodynamics(j,g,i,nts,Data.L,Data.selection, Options);
        if ~isempty(m),
            vol_micro_dyn(:,j,g) = m;
        end
        
        % Computing averaged micro dynamics surface area
        [m,~] = surfacemicrodynamics(j,g,i,nts,Data.L,Data.selection, Options);
        if ~isempty(m),
            surf_micro_dyn(:,j,g) = m;
        end
    end
end


%-------------
% Computing artificial cell lineages
%-------------

% we first construct an array list_cells which will contain the features of
% individual cells - each row correspond to a cell and the column
% correspond to the following cell features: 1) cell population 2)
% generation 3) id mother 4) cell cycle length 5) initial time step 6)
% final time step 7) log daughter/mother volume ratio 8) volume 9) log
% daughter/mother surface ratio 10) surface area

% we need to keep track of the position of the volume and surface area in
% the list of cells
ind_volume = 8;
ind_surface = 10;

% loop over all draws, i.e. all randomly generated cell lineages
for l = 1:ndraw,
    % initialization of the array that will contain all the cells of the
    % current draw
    list_cells = zeros(1,10);
    
    % loop over all the populations
    for j = 1:4,
        
        % loop over the generations
        for k = 0:(ng(j)+1),
            
            g = g0(j)+k;
            Ncells = 2^(g-g0(j))*(n(j));
            
            % initialization of the array that will contain the cells of
            % the current generation and population
            
            % we begin by defining a set of cells without any feature
            list_cells_temp = zeros(Ncells,10);
            
            % their population is defined
            list_cells_temp(:,1) = j;
            % and their generation
            list_cells_temp(:,2) = g;
            
            % we find the indice of the previous group of cells (in terms
            % of generation)
            indmother = find(list_cells(:,1) == j & list_cells(:,2) == g-1);
            % if we can find mothers in sufficient number, i.e. if it is
            % not the initial generation
            if (length(indmother) == 2^(g-g0(j)-1)*(n(j))),
                % we relate the current cells to mothers
                list_cells_temp(:,3) = [indmother;indmother];
                
                % the initial time step is the final time step of their
                % mother
                list_cells_temp(:,5) = list_cells(list_cells_temp(:,3),6);
                
                % we check if the statistics of the cell cycle length have
                % been computed in the embryo, which means that it is not
                % the final generation
                if (pop_parameters.meancyclelength(g,j,1,i)>0),
                    % individual cell features are randomly drawn from the
                    % cell population statistics
                    
                    % cell cycle randomly generated
                    list_cells_temp(:,4) = randn(Ncells,1)*sqrt(pop_parameters.meancyclelength(g,j,2,i))+pop_parameters.meancyclelength(g,j,1,i);
                    
                    % the division time is obtained as the sum of the cell
                    % cycle length and the initial time step of the cell
                    list_cells_temp(:,6) = list_cells_temp(:,4)+list_cells_temp(:,5);
                    
                else
                    % the final time step is the end of the observation
                    % window i.e. tfin
                    list_cells_temp(:,6) = tfin;
                    
                    % the cell cycle length is truncated and defined by the
                    % end of the observation window
                    list_cells_temp(:,4) = list_cells_temp(:,6) - list_cells_temp(:,5);
                    
                end
                % log daughter/mother volume ratio randomly generated
                list_cells_temp(:,7) = pop_parameters.meanlogcoefvol(g,j,1,i)+ sqrt(pop_parameters.meanlogcoefvol(g,j,2,i))*randn(Ncells,1);
                
                % the volume is obtained as the product of the
                % daughter/mother ratio and the volume of the mother
                % cell
                list_cells_temp(:,8) = exp(list_cells_temp(:,7)).*list_cells(list_cells_temp(:,3),8);
                
                % log daughter/mother surface area ratio randomly
                % generated
                list_cells_temp(:,9) = pop_parameters.meanlogcoefsurf(g,j,1,i)+sqrt(pop_parameters.meanlogcoefsurf(g,j,2,i))*randn(Ncells,1);
                
                % the surface area is obtained as the product of the
                % daughter/mother surface ratio and the surface area of
                % the mother cell
                list_cells_temp(:,10) = exp(list_cells_temp(:,9)).*list_cells(list_cells_temp(:,3),10);
                

            else
                % here we handle the case of the initial generation where the
                % beginning of the cell cycles cannot be observed
                % the initial time step is the beginning of the observation
                % window
                list_cells_temp(:,5) = tinit;
                
                % the division times are randomly generated
                list_cells_temp(:,6) = randn(Ncells,1)*sqrt(pop_parameters.meandivisiontime(g,j,2,i))+pop_parameters.meandivisiontime(g,j,1,i);
                
                % the cell cycle length is determined by the difference
                % between the initial and final time steps
                list_cells_temp(:,4) = list_cells_temp(:,6) - list_cells_temp(:,5);
                
                % the volume is randomly generated
                list_cells_temp(:,8) = exp(pop_parameters.meanlogvolume(g,j,1,i) + sqrt(pop_parameters.meanlogvolume(g,j,2,i))*randn(Ncells,1));
                
                % the surface area is randomly generated
                list_cells_temp(:,10) = exp(pop_parameters.meanlogsurface(g,j,1,i) + sqrt(pop_parameters.meanlogsurface(g,j,2,i))*randn(Ncells,1));
            end
            
            % the individual cell features of the current group of cells
            % are concatenated to the previous cells
            list_cells = [list_cells; list_cells_temp];
        end
    end
    
    % now that we have a complete list of cells with realistic individual
    % cell features, we would like to generate the associated cell lineage
    % the cell lineage is stored in a variable called selection_rand which
    % describe the state of each cell at each time step
    
    % this function generate a lineage from individual cell features ...
    selection_rand = selectionfromlistcells(list_cells, ind_volume, ind_surface, time);
    
    % ... however it does not include individual variations undergone by
    % individual cell volume and surface along a cell cycle - we model this
    % in the following couple of loops
    
    % for each population
    for j = 1:4,
        % for each generation
        for g = 6:10,
            % if an averaged cell volume micro dynamics have been defined
            % for the current group of cells (see beginning of this
            % function)
            if sum(abs(vol_micro_dyn(:,j,g)))>0,
                
                % extracting the cells corresponding to the current
                % population and generation
                indtemp1 = find(list_cells(:,1) == j & list_cells(:,2) == g);
                
                % for each cell
                for kl = 1:length(indtemp1),
                    % we find the indices corresponding to the different
                    % stages of the cell in the selection
                    indtemp2 = find(selection_rand(:,1) == indtemp1(kl));
                    % we interpolate the volume micro dynamics on a vector
                    % of the same size as the cell
                    mprime = interp1(1:nts,vol_micro_dyn(:,j,g),(1:((nts-1)/(length(indtemp2)-1)):nts));
                    selection_rand(indtemp2,6) = selection_rand(indtemp2,6).*mprime';
                    % similarly for the surface area
                    mprime = interp1(1:nts,surf_micro_dyn(:,j,g),(1:((nts-1)/(length(indtemp2)-1)):nts));
                    selection_rand(indtemp2,7) = selection_rand(indtemp2,7).*mprime';
                end
                
            end
        end
    end
    
    % computing and keeping track of the embryo level dynamics - Total
    % number of cells, cell volume and cell surface area
    
    for k = 1:length(time),
        for j = 1:4,
            nbcellsstat(k,j,l) = length(find(list_cells(:,1) == j & list_cells(:,5)<=time(k) & list_cells(:,6)>time(k)));
            volcellsstat(k,j,l) = sum(selection_rand(find(selection_rand(:,2) == k & selection_rand(:,3) == j),6));
            surfcellsstat(k,j,l) = sum(selection_rand(find(selection_rand(:,2) == k & selection_rand(:,3) == j),7));
            
        end
        nbcellsstat(k,5,l) = length(find(list_cells(:,5)<=time(k) & list_cells(:,6)>time(k)));
        volcellsstat(k,5,l) = sum(selection_rand(find(selection_rand(:,2) == k),6));
        surfcellsstat(k,5,l) = sum(selection_rand(find(selection_rand(:,2) == k),7));
        
    end
    
end





% Here we compute the average and standard deviation of the generated
% embryo level dynamics for comparison with real data

nbcellsstatmean = zeros(length(time),5,2);
volcellsstatmean = zeros(length(time),5,2);
surfcellsstatmean = zeros(length(time),5,2);

% loop on each time step
for k = 1:length(time),
    % loop on each population and j = 5 correspond to the whole embryo
    for j = 1:5,
        % number of cells
        nbcellsstatmean(k,j,1) = mean(nbcellsstat(k,j,:));
        nbcellsstatmean(k,j,2) = std(nbcellsstat(k,j,:));
        
        % volume
        volcellsstatmean(k,j,1) = mean(volcellsstat(k,j,:));
        volcellsstatmean(k,j,2) = std(volcellsstat(k,j,:));
        
        % surface area
        surfcellsstatmean(k,j,1) = mean(surfcellsstat(k,j,:));
        surfcellsstatmean(k,j,2) = std(surfcellsstat(k,j,:));
        
    end
end


end

