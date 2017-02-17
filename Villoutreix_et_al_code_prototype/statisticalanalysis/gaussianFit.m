function [] = gaussianFit( Data, Options )
%GAUSSIANFIT Computes and draws p-values for each of the distributions of
%individual cell features

% intialization of constants and variables
n = Options.ngen6;
g0 = Options.g0;
lifelengthmin = 20;
I = [1,2,3,4,5];


% pval is an array that contains the value of the p-values for each
% embryo, population, generation, cell feature (division time, cell cycle
% length, volume, surface, log daughter/mother volume, log daughter/mother
% surface)
pval = zeros(5,4,10,6);

divisiontimes=[];
cellcyclelength=[];
volumes=[];
surfaces=[];
daughter_mother_volume=[];
daughter_mother_surface=[];

distributions = cell(5,4,10,6);

% significance threshold
pmin = 0.01;
% lowpval keeps track of the distributions with p<pmin
lowpval = [];

% the first loop goes through the whole cohort of embryos
for i = I,
    % length of a time step
    dt = Data.rescaledtime{i}(2) - Data.rescaledtime{i}(1);
    % the second loop goes through the four cell populations
    for j = 1:4,
        % first generation
        g = g0(i,j);
        k = g-6;
        
        % cells belonging to the current population, generation are
        % extracted, observation of an ending mitosis is required
        indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,6) == 1);
        
        % the considered group of cells should be in a significant number (
        % more than 95%) of the expected number
        if (length(indtemp) >= 0.95*(2^(k)*n(j))),
            
            volumes=[];
            surfaces=[];
            
            % we first consider the distribution of division times
            divisiontimes = (Data.rescaledtime{i}(Data.L{i}(indtemp,4)))/60;
            
            % if the considered cells are observed for a period of more
            % than 20 min, we look at their volume and surface area
            if ((mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60))>lifelengthmin),
                volumes = log(Data.L{i}(indtemp,9));
                surfaces = log(Data.L{i}(indtemp,10));
            end
            
            if ~isempty(divisiontimes),
                l = 1;
                
                % we use Matlab's function chi2gof to compute the p-value
                % for the chi2 goodness of fit test
                [~,p] = chi2gof(divisiontimes);
                
                % we keep track of the division time distribution
                distributions{i,j,g,l} = divisiontimes;
                
                % we keep track of the p-value
                pval(i,j,g,l) = p;
                
                % we keep track of the distribution if the p-value is under
                % the significance threshold
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
            end
            
            if ~isempty(volumes),
                l = 3;
                % computing p-value from chi2 goodness of test
                [~,p] = chi2gof(volumes);
                % keeping track of the distribution
                distributions{i,j,g,l} = volumes;
                % keeping track of the p-value
                pval(i,j,g,l) = p;
                % keeping track of the distribution if the p-value is below
                % pmin
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
            end
            
            if ~isempty(surfaces),
                l = 4;
                % computing p-value from chi2 goodness of test
                [~,p] = chi2gof(surfaces);
                % keeping track of the distribution
                distributions{i,j,g,l} = surfaces;
                % keeping track of the p-value
                pval(i,j,g,l) = p;
                % keeping track of the distribution if the p-value is below
                % pmin
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
                
            end
            
            divisiontimes=[];
            volumes=[];
            surfaces=[];
            
        end
        
        % for subsequent generations
        for g = g0(i,j)+1:10,
            k = g - 6;
            
            % we consider cells of the current generation where the initial
            % and final mitosis are observed
            indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1& Data.L{i}(:,6) == 1);
            
            % if the cells are in significant number (more than 95% of the
            % expected number of cells)
            if (length(indtemp) >= 0.95*(2^(k)*n(j))),
                % distribution of division times
                divisiontimes = (Data.rescaledtime{i}(Data.L{i}(indtemp,4)))/60;
                % distribution of cell cycle lengths
                cellcyclelength = (Data.rescaledtime{i}(Data.L{i}(indtemp,4)+1) - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60;
                % distribution of volumes
                volumes = log(Data.L{i}(indtemp,9));
                % distribution of surfaces
                surfaces = log(Data.L{i}(indtemp,10));
                
                % identifying the mother of the considered cells
                indmother = zeros(length(indtemp),1);
                for m = 1:length(indtemp),
                    indmother(m) =  find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                end
                % if the mother cells are in significant number and with an
                % averaged observation period greater than 20 min
                if (length(indmother) >= 0.95*(2^(k-1)*n(j))) & (mean((Data.rescaledtime{i}(Data.L{i}(indmother,4))+dt - Data.rescaledtime{i}(Data.L{i}(indmother,3)))/60)>lifelengthmin),
                    %log daughter/mother volume coefficient
                    daughter_mother_volume = [];
                    daughter_mother_volume = log(Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                    
                    %log daughter/mother surface coefficient
                    daughter_mother_surface = [];
                    daughter_mother_surface = log(Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                end
                
            else
                % we consider cells of the current population and
                % generation without requiring the observation of the final
                % mitosis
                indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1);
                % if the number of cells is significant and the averaged
                % length of observation is greater than 20 min
                if (length(indtemp) >= 0.95*(2^(k)*n(j)) & (mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60)>lifelengthmin)),
                    % we consider the distribution of volume and surface
                    % area
                    volumes = log(Data.L{i}(indtemp,9));
                    surfaces = log(Data.L{i}(indtemp,10));
                    
                    % identifying the mother of the considered cells
                    indmother = zeros(length(indtemp),1);
                    for m = 1:length(indtemp),
                        indmother(m) =  find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                    end
                    % if the mother cells are in significant number and with an
                    % averaged observation period greater than 20 min
                    if (length(indmother) >= 0.95*(2^(k-1)*n(j))) &(mean((Data.rescaledtime{i}(Data.L{i}(indmother,4))+dt - Data.rescaledtime{i}(Data.L{i}(indmother,3)))/60)>lifelengthmin),
                        %log daughter/mother volume coefficient
                        daughter_mother_volume = log(Data.L{i}(indtemp,9)./Data.L{i}(indmother,9));
                        %log daughter/mother surface coefficient
                        daughter_mother_surface = log(Data.L{i}(indtemp,10)./Data.L{i}(indmother,10));
                    end
                end
            end
            
            % we compute below the p-value for each cell feature - keeping
            % track of the distribution in the structure distributions,
            % p-values in the array pval, and indices of the distribution
            % in the array lowpval if p-value is less than pmin
            if ~isempty(divisiontimes),
                l = 1;
                [~,p] = chi2gof(divisiontimes);
                
                distributions{i,j,g,l} = divisiontimes;
                pval(i,j,g,l) = p;
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
                
            end
            if ~isempty(cellcyclelength),
                l = 2;
                [~,p] = chi2gof(cellcyclelength);
                distributions{i,j,g,l} = cellcyclelength;
                pval(i,j,g,l) = p;
                
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
                
            end
            if ~isempty(volumes),
                l = 3;
                [~,p] = chi2gof(volumes);
                distributions{i,j,g,l} = volumes;
                pval(i,j,g,l) = p;
                
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
                
            end
            if ~isempty(surfaces),
                l = 4;
                [~,p] = chi2gof(surfaces);
                distributions{i,j,g,l} = surfaces;
                pval(i,j,g,l) = p;
                
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
                
            end
            if ~isempty(daughter_mother_volume),
                l = 5;
                [~,p] = chi2gof(daughter_mother_volume);
                distributions{i,j,g,l} = daughter_mother_volume;
                pval(i,j,g,l) = p;
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
            end
            if ~isempty(daughter_mother_surface),
                l = 6;
                [~,p] = chi2gof(daughter_mother_surface);
                distributions{i,j,g,l} = daughter_mother_surface;
                pval(i,j,g,l) = p;
                if p<pmin,
                    lowpval = [lowpval;[i,j,g,l]];
                end
            end
            divisiontimes=[];
            cellcyclelength=[];
            volumes=[];
            surfaces=[];
            daughter_mother_volume=[];
            daughter_mother_surface=[];
        end
        
    end
end


%% Visualization of the distribution of p-values

% we construct an array containing all the p-values different from 0 and
% NaN

X = [];
for i = 1:5,
    for j = 1:4,
        for g = g0(i,j):10,
            for l = 1:6,
                if (pval(i,j,g,l))<1 & ~isnan(pval(i,j,g,l)) & (pval(i,j,g,l))>0,
                    X = [X,pval(i,j,g,l)];
                end
            end
        end
    end
end

% this array is plotted as an histogram with bins of length equal to 0.05
figure(9),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
edges = [0:0.05:1];
[N,BIN] = histc(X,edges);
h = bar(edges,N,'histc');
set(h,'facecolor',[0.25 0.25 0.25]);
set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('p-value');
ylabel('Number of distributions');
title('All statistics');
xlim([0 1])
ylim([0 22])

figure(9); saveas(gcf,'figures/p_values');

%% visualization of the 12 distributions with p-values less than 0.01

type = {'M', 'X', 'logV', 'logS', 'logA', 'logB'};

figure(10),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
c= 1;

% loop on the generations
for g = 7:10,
    
    % loop on the populations
    for j = 1:4,
        ind1 = find(lowpval(:,3) == g & lowpval(:,2) == j);
        
        % loop on the cell features: cell cycle length, division time,
        % volume, surface area, daughter/mother volume/surface area ratio
        for l = [2,1,3,4,5,6],
            
            % extracting the indices of the distributions of the current
            % cell feature with low p-value
            ind2 = find(lowpval(ind1,4) == l);
            
            % loop on these distributions
            for k = 1:length(ind2),

                i = lowpval(ind1(ind2(k)),1);
                j = lowpval(ind1(ind2(k)),2);
                g = lowpval(ind1(ind2(k)),3);
                
                % extracting the distribution
                X = distributions{i,j,g,l};
                % defining the domain
                XX = [min(X):abs(((max(X) - min(X))/10)):max(X)];
                
                % computing the histogram on the domain
                [N,BIN] = histc(X,XX);
                
                % taking care of the last value
                N(end-1) = N(end-1)+N(end);
                N(end) = 0;
                
                % transforming the abscisses values in to double with two
                % decimals and keeping only every other one for
                % visualization purposes
                indicesabs = round(XX(1:2:11)*100)/100;
                
                % plot of the result
                h = subplot(3,4,c);
                bar(XX,N,'histc');
                set(gca, 'FontSize', 14, 'fontName','Arial');
                set(gca,'XTick',indicesabs);
                xlim([min(round(XX*100)/100) max(round(XX*100)/100)])
                title({[type{l},' Emb ',num2str(i), ' ',Options.ident{j}, ' gen ', num2str(g)];[' pval ', num2str(pval(i,j,g,l))]});
                h1 = findobj(gca,'Type','patch');
                set(h1,'FaceColor',[0.5 0.5 0.5],'EdgeColor','w');
                
                
                c = c+1;
            end
        end
    end
end

figure(10); saveas(gcf,'figures/worst_fits');

%% visualization of 12 distributions with the highest p-value

type = {'M', 'X', 'logV', 'logS', 'logA', 'logB'};

figure(11),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
c = 1;

% the first loop is on the cell features
for l = [2,1,3,4,5,6],
    % we consider the p-values corresponding to the current cell feature
    temp = pval(:,:,:,l);
    % we draw the distributions with the two highest p-values
    for k = 1:2,
        % find the distribution with the highest p-value
        [C, indi] = max(temp(:));
        [i,j,g] = ind2sub(size(temp),indi);
        
        % extracting the distribution
        X = distributions{i,j,g,l};
        
        % defining the domain
        XX = [min(X):abs(((max(X) - min(X))/10)):max(X)];
        
        % computing the histogram on the domain
        [N,BIN] = histc(X,XX);
        
        % taking care of the last value
        N(end-1) = N(end-1)+N(end);
        N(end) = 0;
        
        % transforming the abscisses values in to double with two
        % decimals and keeping only every other one for visualisation
        % purposes
        indicesabs = round(XX(1:2:11)*100)/100;
        
        % plot of the result
        h = subplot(3,4,c);
        bar(XX,N,'histc')
        set(gca, 'FontSize', 14, 'fontName','Arial');
        set(gca,'XTick',indicesabs);
        xlim([min(round(XX*100)/100) max(round(XX*100)/100)])
        title({[type{l},' Emb ',num2str(i), ' ',Options.ident{j}, ' gen ', num2str(g)];[' pval ', num2str(pval(i,j,g,l))]});
        h1 = findobj(gca,'type','patch');
        set(h1,'FaceColor',[0.3 0.3 0.3],'EdgeColor','w');
        
        % removing the current distribution from the list of p-value for
        % the second loop
        temp(i,j,g) = 0;
        
        c = c+1;
    end
end

figure(11); saveas(gcf,'figures/best_fits');

end

