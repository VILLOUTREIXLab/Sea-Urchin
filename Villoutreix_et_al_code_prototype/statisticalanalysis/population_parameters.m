function [ pop_parameters ] = population_parameters(Data, Options)
%POPULATION_PARAMETERS computes the statistical parameters for each cell
% group and each individual cell feature: division time, cell cycle length,
% volume, surface area, mother/daughter volume/surface area coefficients
 
% Initialization of the variables
Ng = 12;
Nsp = 4;
meancyclelength = zeros(Ng, Nsp,2,5);
meandivisiontime = zeros(Ng, Nsp, 2,5);
meanvolume = zeros(Ng, Nsp,2,5);
meancoefvol = zeros(Ng, Nsp, 2,5);
meanlogvolume = zeros(Ng, Nsp,2,5);
meanlogcoefvol = zeros(Ng, Nsp, 2,5);
meanlogsurface = zeros(Ng, Nsp,2,5);
meanlogcoefsurf = zeros(Ng, Nsp, 2,5);
meansurface = zeros(Ng, Nsp,2,5);
meancoefsurf = zeros(Ng, Nsp,2,5);

% number of cells for each population at generation 6
n = Options.ngen6;
% initial generation number for each population of each embryo
g0 = Options.g0;

% Minimal observed cell cycle length for a cell to be considered significant
cyclelengthmin = 20;

% loop over all the embryos of the cohort
for i = 1:Options.Nb,
    
    % legnth of a time step
    dt = Data.rescaledtime{i}(2) - Data.rescaledtime{i}(1);
    
    % loop over all the cell populations
    for j = 1:4,
        % for the initial generation
        g = g0(i,j);
        k = g-6;
        % we extract the cells belonging to the considered population,
        % generation, given embryo, and with observed mitosis at the end of
        % the cycle
        indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,6) == 1);
        % if the number of extracted cells is considered significant (95%
        % of the expected number of cells)
        if (length(indtemp) >= 0.95*(2^(k)*n(j)))
            %divisiontime
            X = [];
            X = Data.rescaledtime{i}(Data.L{i}(indtemp,4))/60;
            meandivisiontime(g,j,1,i) = mean(X);
            meandivisiontime(g,j,2,i) = var(X);
            
            % if the averaged length of  cell cycles in the set of
            % extracted cells are greater than 20 min
            if  ((mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60))>cyclelengthmin),
                % computing mean and variance for ...
                % ... log volume
                X = [];
                X = log(Data.L{i}(indtemp,9));
                meanlogvolume(g,j,1,i) = mean(X);
                meanlogvolume(g,j,2,i) = var(X);
                
                % ... volume
                X = [];
                X = (Data.L{i}(indtemp,9));
                meanvolume(g,j,1,i) = mean(X);
                meanvolume(g,j,2,i) = var(X);
                
                % ... log surface area
                X = [];
                X = log(Data.L{i}(indtemp,10));
                meanlogsurface(g,j,1,i) = mean(X);
                meanlogsurface(g,j,2,i) = var(X);
                
                % ... surface area
                X = [];
                X = (Data.L{i}(indtemp,10));
                meansurface(g,j,1,i) = mean(X);
                meansurface(g,j,2,i) = var(X);
            end
        end
        
        % for subsequent generations
        for g = (g0(i,j)+1):10,
            k = g - 6;
            
            % we extract cells belonging to the considered population,
            % generation, given embryo, and with observed initial and final
            % division
            indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1 & Data.L{i}(:,6) == 1);
            % if the number of extracted cells is considered significant (95%
            % of the expected number of cells)
            
            if (length(indtemp) >= 0.95*(2^(k)*n(j))),
                % computing mean and variance for ...
                % ... division time
                X = [];
                X = Data.rescaledtime{i}(Data.L{i}(indtemp,4))/60;
                meandivisiontime(g,j,1,i) = mean(X);
                meandivisiontime(g,j,2,i) = var(X);
                
                % ... life length
                X = [];
                X = (Data.rescaledtime{i}(Data.L{i}(indtemp,4)+1) - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60;
                meancyclelength(g,j,1,i) = mean(X);
                meancyclelength(g,j,2,i) = var(X);
                
                % ... log volume
                X = [];
                X = log(Data.L{i}(indtemp,9));
                meanlogvolume(g,j,1,i) = mean(X);
                meanlogvolume(g,j,2,i) = var(X);
                
                % ... volume
                X = [];
                X = (Data.L{i}(indtemp,9));
                meanvolume(g,j,1,i) = mean(X);
                meanvolume(g,j,2,i) = var(X);
                
                % ... log surface area
                X = [];
                X = log(Data.L{i}(indtemp,10));
                meanlogsurface(g,j,1,i) = mean(X);
                meanlogsurface(g,j,2,i) = var(X);
                
                % ... surface area
                X = [];
                X = (Data.L{i}(indtemp,10));
                meansurface(g,j,1,i) = mean(X);
                meansurface(g,j,2,i) = var(X);
                
                % finding the indices of the mother for each cells
                indmother = zeros(length(indtemp),1);
                for m = 1:length(indtemp),
                    indmother(m) =  find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                end
                
                % if the set of mother cells is in significant numer (95%
                % criteria) and with significant length (at least 20 min)
                if (length(indmother) >= 0.95*(2^(k-1)*n(j))) & (mean((Data.rescaledtime{i}(Data.L{i}(indmother,4))+dt - Data.rescaledtime{i}(Data.L{i}(indmother,3)))/60)>cyclelengthmin),
                    % computing mean and variance for ...
                    % ... log daughter/mother volume ratio
                    X = [];
                    X = log(Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                    if (~isnan(X)),
                        meanlogcoefvol(g,j,1,i) = mean(X);
                        meanlogcoefvol(g,j,2,i) = var(X);
                    end
                    
                    % ... daughter/mother volume ratio
                    X = [];
                    X = (Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                    if (~isnan(X)),
                        meancoefvol(g,j,1,i) = mean(X);
                        meancoefvol(g,j,2,i) = var(X);
                    end
                    
                    % ... log daughter/mother surface ratio
                    X = [];
                    X = log(Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                    if (~isnan(X)),
                        meanlogcoefsurf(g,j,1,i) = mean(X);
                        meanlogcoefsurf(g,j,2,i) = var(X);
                    end
                    % ... daughter/mother surface ratio
                    X = [];
                    X = (Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                    if (~isnan(X)),
                        meancoefsurf(g,j,1,i) = mean(X);
                        meancoefsurf(g,j,2,i) = var(X);
                    end
                    
                end
            else
                % looking at incomplete cell cycles, where only the initial
                % mitosis can be observed (last cycle of the observation
                % window
                indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1);
                % if the considered set of cells validate the significance
                % requirements (number of cells, and length of cell cycle)
                if (length(indtemp) >= 0.95*(2^(k)*n(j)) & (mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60)>cyclelengthmin)),
                    % computing mean and variance for ...
                    % ... logvolume
                    X = [];
                    X = log(Data.L{i}(indtemp,9));
                    meanlogvolume(g,j,1,i) = mean(X);
                    meanlogvolume(g,j,2,i) = var(X);
                    % ... volume
                    X = [];
                    X = (Data.L{i}(indtemp,9));
                    meanvolume(g,j,1,i) = mean(X);
                    meanvolume(g,j,2,i) = var(X);
                    % ... log surface
                    X = [];
                    X = log(Data.L{i}(indtemp,10));
                    meanlogsurface(g,j,1,i) = mean(X);
                    meanlogsurface(g,j,2,i) = var(X);
                    % ... surface
                    X = [];
                    X = (Data.L{i}(indtemp,10));
                    meansurface(g,j,1,i) = mean(X);
                    meansurface(g,j,2,i) = var(X);
                    % finding the indices of the mother for each cells
                    indmother = zeros(length(indtemp),1);
                    for m = 1:length(indtemp),
                        indmother(m) =  find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                    end
                    % if the set of mother cells is in significant numer (95%
                    % criteria) and with significant length (at least 20 min)
                    if (length(indmother) >= 0.95*(2^(k-1)*n(j))) & (mean((Data.rescaledtime{i}(Data.L{i}(indmother,4))+dt - Data.rescaledtime{i}(Data.L{i}(indmother,3)))/60)>cyclelengthmin),
                        % computing mean and variance for ...
                        % ... log daughter/mother volume ratio
                        X = [];
                        X = log(Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                        if (~isnan(X)),
                            meanlogcoefvol(g,j,1,i) = mean(X);
                            meanlogcoefvol(g,j,2,i) = var(X);
                        end
                        % ... daughter/mother volume ratio
                        X = [];
                        X = (Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                        if (~isnan(X)),
                            meancoefvol(g,j,1,i) = mean(X);
                            meancoefvol(g,j,2,i) = var(X);
                        end
                        % ... log daughter/mother surface ratio
                        X = [];
                        X = log(Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                        if (~isnan(X)),
                            meanlogcoefsurf(g,j,1,i) = mean(X);
                            meanlogcoefsurf(g,j,2,i) = var(X);
                        end
                        % ... daughter/mother surface ratio
                        X = [];
                        X = (Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                        if (~isnan(X)),
                            meancoefsurf(g,j,1,i) = mean(X);
                            meancoefsurf(g,j,2,i) = var(X);
                        end
                    end
                end
            end
        end
    end
end

% values of the statistics of cell features are stored in the object
% pop_parameters which is returned by the function

pop_parameters.meancyclelength = meancyclelength;
pop_parameters.meandivisiontime = meandivisiontime;
pop_parameters.meanvolume = meanvolume;
pop_parameters.meancoefvol = meancoefvol;
pop_parameters.meanlogvolume = meanlogvolume;
pop_parameters.meanlogcoefvol = meanlogcoefvol;
pop_parameters.meanlogsurface = meanlogsurface;
pop_parameters.meanlogcoefsurf = meanlogcoefsurf;
pop_parameters.meansurface = meansurface;
pop_parameters.meancoefsurf = meancoefsurf;

end

