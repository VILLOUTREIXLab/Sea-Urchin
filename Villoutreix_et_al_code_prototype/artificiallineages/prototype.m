function [] = prototype(ndraw, Data, Options, proto_parameters)
%PROTOTYPE generate random cell lineages based on prototypical
%statistics and draws embryo-level dynamics based on them


% Initial stage - 32 cells
n = Options.ngen6;
g0 = [6,6,6,6];
% number of cycles for each population
ng = [0,2,3,3];

% computing a time vector large enough
tmini_ = [];
tmax_ = [];
for i = 1:5,
    tmini_(i) = Data.rescaledtime{i}(1)/60;
    tmax_(i) = Data.rescaledtime{i}(Options.tmax(i))/60;
end
tmini = min(tmini_);
tmaxi = max(tmax_);
time = (tmini:tmaxi);

% initialization of the arrays that will keep track of the embryo-level
% dynamics from the random cell lineages
nbcellsstat = zeros(length(time),5,ndraw);
volcellsstat = zeros(length(time),5,ndraw);
surfcellsstat = zeros(length(time),5,ndraw);

%-------------
% Computing volume and surface area micro dynamics as an averaged over the
% cohort
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
        X = [];
        % for the volume
        % loop on the embryos
        for i = 1:Options.Nb,
        % Computing averaged micro dynamics volume for embryo i
        [m,~] = volumemicrodynamics(j,g,i,nts,Data.L,Data.selection, Options);
        % stored in an array
        X = [X,m];
        end
        % the average over the cohort is computed and stored
        if ~isempty(X),
           vol_micro_dyn(:,j,g) = mean(X,2); 
        end
        
        % for the surface area
        % loop on the embryos
        for i = 1:Options.Nb,
        % Computing averaged micro dynamics volume for embryo i
        [m,~] = surfacemicrodynamics(j,g,i,nts,Data.L,Data.selection, Options);
        % stored in an array
        X = [X,m];
        end
        % the average over the cohort is computed and stored
        if ~isempty(X),
           surf_micro_dyn(:,j,g) = mean(X,2); 
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

    % 1) Sous Pop 2) Generation 3) Num Mother 4) Duree de vie 5) moment debut
    % 6) moment fin 7) volume 8) surface

% we need to keep track of the position of the volume and surface area in
% the list of cells
ind_volume = 7;
ind_surface = 8;

% loop over all draws, i.e. all randomly generated cell lineages
for l = 1:ndraw,
    % initialization of the array that will contain all the cells of the
    % current draw
    listecells = zeros(1,9);

    
    % loop over the populations
    for j = 1:4,
        % loop over the generations
        for k = 0:(ng(j)+1),
            
            g = g0(j)+k;
            Ncells = 2^(g-g0(j))*(n(j));
            
            % initialization of the array that will contain the cells of
            % the current generation and population
            
            % we begin by defining a set of cells without any feature
            list_cells_temp = zeros(Ncells,9);
            
            % their population is defined
            list_cells_temp(:,1) = j;
            % and their generation
            list_cells_temp(:,2) = g;
            
            % we find the indice of the previous group of cells (in terms
            % of generation)
            indmother = find(listecells(:,1) == j & listecells(:,2) == g-1);
            % if we can find mothers in sufficient number, i.e. if it is
            % not the initial generation
            if (length(indmother) == 2^(g-g0(j)-1)*(n(j))),
                % we relate the current cells to mothers
                list_cells_temp(:,3) = [indmother;indmother];
                
                % the initial time step is the final time step of their
                % mother
                list_cells_temp(:,5) = listecells(list_cells_temp(:,3),6);
                
                % we check if the statistics of the cell cycle length have
                % been computed for the prototype, which means that it is 
                % not the final generation
                if proto_parameters.meancyclelengthProto(g,j,1)>0,
                    % individual cell features are randomly drawn from the
                    % prototypicam cell population statistics
                    
                    % cell cycle randomly generated
                    list_cells_temp(:,4) = randn(Ncells,1)*sqrt(proto_parameters.meancyclelengthProto(g,j,2))+proto_parameters.meancyclelengthProto(g,j,1);
                    
                    % the division time is obtained as the sum of the cell 
                    % cycle length and the initial time step of the cell
                    list_cells_temp(:,6) = list_cells_temp(:,4)+list_cells_temp(:,5);
                    
                else
                    % the final time step is the end of the observation
                    % window i.e. tmaxi
                    list_cells_temp(:,6) = tmaxi;
                    
                    % the cell cycle length is truncated and defined by the
                    % end of the observation window
                    list_cells_temp(:,4) = list_cells_temp(:,6) - list_cells_temp(:,5);
                end
                % the volume is obtained as the product of the randomly
                % generated daughter/mother ratio and the volume of the 
                % mother cell
                list_cells_temp(:,7) = listecells(list_cells_temp(:,3),7).*exp(randn(Ncells,1)*sqrt(proto_parameters.meanlogcoefvolProto(g,j,2))+proto_parameters.meanlogcoefvolProto(g,j,1));
                
                % the surface area is obtained as the product of the randomly
                % generated daughter/mother ratio and the volume of the 
                % mother cell
                list_cells_temp(:,8) = listecells(list_cells_temp(:,3),8).*exp(randn(Ncells,1)*sqrt(proto_parameters.meanlogcoefsurfProto(g,j,2))+proto_parameters.meanlogcoefsurfProto(g,j,1));
                
                
            else
                % here we handle the case of the initial generation where the
                % beginning of the cell cycles cannot be observed
                % the initial time step is the beginning of the observation
                % window
                list_cells_temp(:,5) = tmini;
                
                % the division times are randomly generated
                list_cells_temp(:,6) = randn(Ncells,1)*sqrt(proto_parameters.meandivisiontimeProto(g,j,2))+proto_parameters.meandivisiontimeProto(g,j,1);
                
                % the cell cycle length is determined by the difference
                % between the initial and final time steps
                list_cells_temp(:,4) = list_cells_temp(:,6) - list_cells_temp(:,5);
                
                % the volume is randomly generated
                list_cells_temp(:,7) = exp(randn(Ncells,1)*sqrt(proto_parameters.meanlogvolumeProto(g,j,2)) + proto_parameters.meanlogvolumeProto(g,j,1));
                
                % the surface area is randomly generated
                list_cells_temp(:,8) = exp(randn(Ncells,1)*sqrt(proto_parameters.meanlogsurfaceProto(g,j,2)) + proto_parameters.meanlogsurfaceProto(g,j,1));
                
            end
            
            % the individual cell features of the current group of cells
            % are concatenated to the previous cells
            listecells = [listecells; list_cells_temp];
        end
    end
    
    % now that we have a complete list of cells with realistic individual
    % cell features, we would like to generate the associated cell lineage
    % the cell lineage is stored in a variable called selection_rand which
    % describe the state of each cell at each time step
    
    % this function generate a lineage from individual cell features ...
    selection_rand = selectionfromlistcells(listecells, ind_volume,ind_surface, time);
    
    
    % ... however it does not include individual variations undergone by
    % individual cell volume and surface along a cell cycle - we model this
    % in the following couple of loops
    
    % for each population
    for j = 1:4,
        % for each generation
        for g = 6:10,
            if sum(abs(vol_micro_dyn(:,j,g)))>0,
                
                % extracting the cells corresponding to the current
                % population and generation
                indtemp1 = find(listecells(:,1) == j & listecells(:,2) == g);
                
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
            nbcellsstat(k,j,l) = length(find(listecells(:,1) == j & listecells(:,5)<=time(k) & listecells(:,6)>time(k)));
            volcellsstat(k,j,l) = sum(selection_rand(find(selection_rand(:,2) == k & selection_rand(:,3) == j),6));
            surfcellsstat(k,j,l) = sum(selection_rand(find(selection_rand(:,2) == k & selection_rand(:,3) == j),7));
        end
        nbcellsstat(k,5,l) = length(find(listecells(:,5)<=time(k) & listecells(:,6)>time(k)));
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




%-------------------
% Visualization
%-------------------

% the three figures, Number of cells, Cell Volume, Cell Surface Area, are 
% generated with the same layout - Comments are provided for the number of
% cells

% defining the domain used to draw gray area for standard deviation
X = [time';flipud(time')];

% Number of cells

figure(22),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% for each generation and for the whole embryo
for j = 1:5,
    % drawing the summary of the randomly generated cell lineages using
    % prototypical statistics
    
    % we draw a gray area corresponding to one standard deviation above the
    % mean number of cells ...
    Y = [nbcellsstatmean(:,j,1);flipud(nbcellsstatmean(:,j,1)+nbcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    hold off,
    % ... and below
    Y = [nbcellsstatmean(:,j,1);flipud(nbcellsstatmean(:,j,1)-nbcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    hold off,
    
    hold on,
    % we draw the mean number of cells
    plot(time/60,nbcellsstatmean(:,j,1),'-','LineWidth',5,'Color',[0.2,0.2,0.2]);
    hold off,
    
    % for each embryo of the cohort we draw the total number of cells
    for i = 1:Options.Nb,
        hold on,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600, Data.NbCells(1:Options.tmax(i),j,i),'-', 'Color',Options.colors_embryos(i,:),'LineWidth',2.5);
        hold off,
    end
    
end

abslim = [4.2 9.5];
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Cell number');
title('Cell number for the prototypical embryo');


figure(22); saveas(gcf,'figures/nbcells_proto');

% Volume
figure(23),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,

    Y = [volcellsstatmean(:,j,1);flipud(volcellsstatmean(:,j,1)+volcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    
    hold off,
    Y = [volcellsstatmean(:,j,1);flipud(volcellsstatmean(:,j,1)-volcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    
    hold off,
    hold on,
    
    plot(time/60,volcellsstatmean(:,j,1),'-','LineWidth',5,'Color',[0.2,0.2,0.2]);
    
    for i = 1:Options.Nb,
        hold on,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600, Data.VolCells(1:Options.tmax(i),j,i),'-', 'Color',Options.colors_embryos(i,:),'LineWidth',2.5);
        hold off,
    end
end

abslim = [4.2 9.5];

set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Volume');
title('Volume of the cells for the prototypical embryo');

figure(23); saveas(gcf,'figures/cellvolume_proto');

% Surface area
figure(24),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
    Y = [surfcellsstatmean(:,j,1);flipud(surfcellsstatmean(:,j,1)+surfcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    hold off,
    Y = [surfcellsstatmean(:,j,1);flipud(surfcellsstatmean(:,j,1)-surfcellsstatmean(:,j,2))];
    hold on,
    fill(X/60,Y,[0.8 0.8 0.8],'EdgeColor','none');
    hold off,
    hold on,
    plot(time/60,surfcellsstatmean(:,j,1),'-','LineWidth',5,'Color',[0.2,0.2,0.2]);
    hold off,
    for i = 1:Options.Nb,
        hold on,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600, Data.SurfCells(1:Options.tmax(i),j,i),'-', 'Color',Options.colors_embryos(i,:),'LineWidth',2.5);
        hold off,
    end
end
xlim([4.2 9.5]);
ylim([0 260000]);

set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim);
ylabel('Surface');
title('Surface of the cells with prototypical embryo');

figure(24); saveas(gcf,'figures/cellsurface_proto');

end

