function [] = correlations( Data, Options )
%CORRELATIONS computes correlation coefficients for distribution of cell
%features in the cell lineage

% Initialization of the constants and variables
I = 1:Options.Nb;
n = Options.ngen6;
g0 = Options.g0;

coef = zeros(10,4,5,6);
pvals = zeros(10,4,5,6);

cyclelengthmin = 20;

% first loop on the embryos
for i = I,
    % duration of a time step
    dt = Data.rescaledtime{i}(2) - Data.rescaledtime{i}(1);
    
    % second loop on the populations
    for j = 1:4,
        
        % for the initial generation
        g = g0(i,j);
        
        %-----------
        % cell cycle length
        %-----------
        
        % we look at the cells with complete cell cycles
        [indtemp] = find(Data.L{i}(:,2) == j & Data.L{i}(:,5) == 1 & Data.L{i}(:,6) == 1 & Data.L{i}(:,7) == g);
        % and their mother cells
        mothers = unique(Data.L{i}(indtemp,8));
        
        indtemp1 = [];
        indtemp2 = [];
        indmother = [];
        
        % loop on the mothers
        for l = 1:length(mothers),
            % find the two daughter cells if they exist
            itemp = indtemp(find(Data.L{i}(indtemp,8) == mothers(l)));
            
            if length(itemp)>=2,
                % cycle length for the two daughter cells
                ll1 = ((Data.rescaledtime{i}(Data.L{i}(itemp(1),4))+dt - Data.rescaledtime{i}(Data.L{i}(itemp(1),3)))/60);
                ll2 = ((Data.rescaledtime{i}(Data.L{i}(itemp(2),4))+dt - Data.rescaledtime{i}(Data.L{i}(itemp(2),3)))/60);
                
                % they are ordered in a common vector
                if ll1>ll2,
                    indtemp1 = [indtemp1;(itemp(1))];
                    indtemp2 = [indtemp2;(itemp(2))];
                else
                    indtemp1 = [indtemp1;(itemp(2))];
                    indtemp2 = [indtemp2;(itemp(1))];
                end
                indmother = [indmother;Data.L{i}(mothers(l),1)];
            end
        end
        
        % if the number of cells having an identified mother is significant
        % (at least 95 % of the expected number of cells)
        if (((length(indtemp1)+length(indtemp2)) >= 0.95*2^(g-6)*n(j))),
            % computing the whole list of daughter cells cycles
            ll1 = ((Data.rescaledtime{i}(Data.L{i}(indtemp1,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp1,3)))/60);
            ll2 = ((Data.rescaledtime{i}(Data.L{i}(indtemp2,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp2,3)))/60);
            
            % computing the correlation between the cell cycle lengths of
            % the pairs of sisters
            [coef_,pvals_] = corrcoef(ll1, ll2);
            coef(g,j,i,1) = coef_(1,2);
            pvals(g,j,i,1) = pvals_(1,2);
        end
        
        % for the subsequent generations
        for g = g0(i,j)+1:10,
            %-----------
            % cell cycle length
            %-----------
            
            % identifying cells of the current population generation with
            % observed initial and final mitoses
            [indtemp] = find(Data.L{i}(:,2) == j & Data.L{i}(:,5) == 1 & Data.L{i}(:,6) == 1 & Data.L{i}(:,7) == g);
            
            %extracting the set of mother cells
            mothers = unique(Data.L{i}(indtemp,8));
            
            indtemp1 = [];
            indtemp2 = [];
            indmother = [];
            
            % loop on the mothers
            for l = 1:length(mothers),
                % find the two daughter cells if they exist
                itemp = indtemp(find(Data.L{i}(indtemp,8) == mothers(l)));
                
                if length(itemp)>=2,
                    % cycle length of the two daughter cells
                    ll1 = ((Data.rescaledtime{i}(Data.L{i}(itemp(1),4))+dt - Data.rescaledtime{i}(Data.L{i}(itemp(1),3)))/60);
                    ll2 = ((Data.rescaledtime{i}(Data.L{i}(itemp(2),4))+dt - Data.rescaledtime{i}(Data.L{i}(itemp(2),3)))/60);
                    
                    % they are ordered in a common vector
                    if ll1>ll2,
                        indtemp1 = [indtemp1;(itemp(1))];
                        indtemp2 = [indtemp2;(itemp(2))];
                    else
                        indtemp1 = [indtemp1;(itemp(2))];
                        indtemp2 = [indtemp2;(itemp(1))];
                    end
                    indmother = [indmother;Data.L{i}(mothers(l),1)];
                end
            end
            
            % if the number of cells having an identified mother is significant
            % (at least 95 % of the expected number of cells)
            if (((length(indtemp)) >= 0.95*2^(g-6)*n(j))),
                % computing the whole list of daughter cells cycles
                ll1 = ((Data.rescaledtime{i}(Data.L{i}(indtemp1,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp1,3)))/60);
                ll2 = ((Data.rescaledtime{i}(Data.L{i}(indtemp2,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp2,3)))/60);
                
                % computing the correlation between the cell cycle lengths of
                % the pairs of sisters
                [coef_,pvals_] = corrcoef(ll1, ll2);
                coef(g,j,i,1) = coef_(1,2);
                pvals(g,j,i,1) = pvals_(1,2);
            end
            
            %-----------
            % daughter/mother volume ratio
            %-----------
            
            % extracting cells belonging to the current population and
            % generation and with observed initial mitosis
            [indtemp] = find(Data.L{i}(:,2) == j & Data.L{i}(:,5) == 1 & Data.L{i}(:,7) == g);
            % extracting the set of mother cells
            mothers = unique(Data.L{i}(indtemp,8));
            
            indtemp1 = [];
            indtemp2 = [];
            indmother = [];
            
            % if the set of cells and the set of mother cells present
            % averaged cell cycle length greater than 20 min
            if ((mean((Data.rescaledtime{i}(Data.L{i}(mothers,4))+dt - Data.rescaledtime{i}(Data.L{i}(mothers,3)))/60) >cyclelengthmin) & (mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60) >cyclelengthmin)),
                % loop on the mothers
                for l = 1:length(mothers),
                    % identify the two daughter cells if they exist
                    itemp = indtemp(find(Data.L{i}(indtemp,8) == mothers(l)));
                    if length(itemp)>=2,
                        % computing daughter/mother volume ratio for the
                        % two daughters
                        ll1 = (Data.L{i}(itemp(1),9)./Data.L{i}(mothers(l),9))';
                        ll2 = (Data.L{i}(itemp(2),9)./Data.L{i}(mothers(l),9))';
                        
                        % ordering them
                        if ll1>ll2,
                            indtemp1 = [indtemp1;(itemp(1))];
                            indtemp2 = [indtemp2;(itemp(2))];
                        else
                            indtemp1 = [indtemp1;(itemp(2))];
                            indtemp2 = [indtemp2;(itemp(1))];
                        end
                        indmother = [indmother;Data.L{i}(mothers(l),1)];
                    end
                end
            end
            
            % if the cells are in significant number (at least 95% of the
            % expected number of cells)
            if (((length(indtemp1)+length(indtemp2)) >= 0.95*2^(g-6)*n(j))),
                % computing the whole list of daughter/mother volume ratio
                ll1 = (Data.L{i}(indtemp1,9)./Data.L{i}(indmother,9))';
                ll2 = (Data.L{i}(indtemp2,9)./Data.L{i}(indmother,9))';
                
                %computing the correlation between the pairs of sisters
                [coef_,pvals_] = corrcoef(ll1, ll2);
                coef(g,j,i,2) = coef_(1,2);
                pvals(g,j,i,2) = pvals_(1,2);
                
            end
            
            %-----------
            % daughter/mother surface area ratio
            %-----------
            
            % extracting cells belonging to the current population and
            % generation and with observed initial mitosis
            [indtemp] = find(Data.L{i}(:,2) == j & Data.L{i}(:,5) == 1 & Data.L{i}(:,7) == g);
            % extracting the set of mother cells
            mothers = unique(Data.L{i}(indtemp,8));
            
            indtemp1 = [];
            indtemp2 = [];
            indmother = [];
            
            % if the set of cells and the set of mother cells present
            % averaged cell cycle length greater than 20 min
            if ((mean((Data.rescaledtime{i}(Data.L{i}(mothers,4))+dt - Data.rescaledtime{i}(Data.L{i}(mothers,3)))/60) >cyclelengthmin) & (mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60) >cyclelengthmin)), 
                % loop on the mothers
                for l = 1:length(mothers),
                    % identify the two daughter cells if they exist
                    itemp = indtemp(find(Data.L{i}(indtemp,8) == mothers(l)));
                    if length(itemp)>=2,
                        % computing daughter/mother surface area ratio for
                        % the two daughters
                        ll1 = (Data.L{i}(itemp(1),10)./Data.L{i}(mothers(l),10))';
                        ll2 = (Data.L{i}(itemp(2),10)./Data.L{i}(mothers(l),10))';
                        
                        % ordering them
                        if ll1>ll2,
                            indtemp1 = [indtemp1;(itemp(1))];
                            indtemp2 = [indtemp2;(itemp(2))];
                        else
                            indtemp1 = [indtemp1;(itemp(2))];
                            indtemp2 = [indtemp2;(itemp(1))];
                        end
                        indmother = [indmother;Data.L{i}(mothers(l),1)];
                    end
                end
            end
            
            % if the cells are in significant number (at least 95% of the
            % expected number of cells)
            if (((length(indtemp1)+length(indtemp2)) >= 0.95*2^(g-6)*n(j))),
                % computing the whole list of daughter/mother surface area
                % ratio
                ll1 = (Data.L{i}(indtemp1,10)./Data.L{i}(indmother,10))';
                ll2 = (Data.L{i}(indtemp2,10)./Data.L{i}(indmother,10))';
                
                %computing the correlation between the pairs of sisters
                [coef_,pvals_] = corrcoef(ll1, ll2);
                coef(g,j,i,3) = coef_(1,2);
                pvals(g,j,i,3) = pvals_(1,2);
            end
            
            %-----------
            %-----------
            % mother/daughter relationships
            %-----------
            %-----------
            
            %-----------
            % daughter's cell cycle length and mother's division time
            %-----------
            
            % extracting cells with observed initial and final mitosis
            indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1 & Data.L{i}(:,6) == 1);
            
            % if they are in significant number (95% rule)
            if (length(indtemp) >= 0.95*(2^(g-6)*n(j))),
                %identifying the position of the mothers in L for each
                %daughter cell
                indmother = zeros(length(indtemp),1);
                for m = 1:length(indtemp),
                    indm = find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                    if ((Data.L{i}(indm,6) == 1)),
                        indmother(m) = indm;
                    end
                end
                
                % if the set of mother is in significant number
                if (length(unique(indmother)) >= 0.95*(2^(g-6-1)*n(j))),
                    % computing the cell cycle length of the daughter
                    X1 = Data.rescaledtime{i}(Data.L{i}(indtemp(find(indmother)),4)+1)-Data.rescaledtime{i}(Data.L{i}(indtemp(find(indmother)),3));
                    % and the division time of the mother
                    X2 = Data.rescaledtime{i}(Data.L{i}(indmother(find(indmother)),4));
                    
                    % computing the correlation coefficient between
                    % daughter and mother features
                    [coef_,pvals_] = corrcoef(X1,X2);
                    coef(g,j,i,4) = coef_(1,2);
                    pvals(g,j,i,4) = pvals_(1,2);
                    
                end
            end
            
            %-----------
            % daughter's volume/surface ratio and mother's volume/surface
            %-----------
            
            % extracting the cells belonging to current population and
            % generation
            indtemp = find(Data.L{i}(:,2) == j & Data.L{i}(:,7) == g & Data.L{i}(:,5) == 1);
            % if cells are in significant number and have averaged cell
            % cycle length greater than 20 min
            if (length(indtemp) >= 0.95*(2^(g-6)*n(j)) & (mean((Data.rescaledtime{i}(Data.L{i}(indtemp,4))+dt - Data.rescaledtime{i}(Data.L{i}(indtemp,3)))/60)>cyclelengthmin)),
                
                % identifying mother cells
                indmother = zeros(length(indtemp),1);
                for m = 1:length(indtemp),
                    indmother(m) =  find(Data.L{i}(:,1) == Data.L{i}(indtemp(m),8));
                end
                
                % if the mothers are in significan number and have averaged
                % cell cycle length greater than 20 min
                if (length(find(unique(indmother))) >= 0.95*(2^(g-6-1)*n(j))) & (mean((Data.rescaledtime{i}(Data.L{i}(indmother,4))+dt - Data.rescaledtime{i}(Data.L{i}(indmother,3)))/60)>cyclelengthmin ),
                    
                    % daughter's daughter/mother volume ratio
                    X1 =  Data.L{i}(indtemp(find(indmother)),9)./Data.L{i}(indmother(find(indmother)),9);
                    % and mother's volume
                    X2 = Data.L{i}(indmother(find(indmother)),9);
                    
                    % computing the correlation
                    [coef_,pvals_] = corrcoef(X1,X2);
                    coef(g,j,i,5) = coef_(1,2);
                    pvals(g,j,i,5) = pvals_(1,2);
                    
                    % daughter's daughter/mother surface area ratio
                    X1 =  Data.L{i}(indtemp(find(indmother)),10)./Data.L{i}(indmother(find(indmother)),10);
                    % and mother's surface area
                    X2 = Data.L{i}(indmother(find(indmother)),10);
                    
                    % computing the correlation
                    [coef_,pvals_] = corrcoef(X1,X2);
                    coef(g,j,i,6) = coef_(1,2);
                    pvals(g,j,i,6) = pvals_(1,2);
                end
            end
        end
    end
end


%% Visualization of the squared correlation coefficients

J = [1,2,3,4];
figure(12),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
X = [];

% loop on the embryos
for i = I,
    % loop on the populations
    for j = J,
        % loop on the generations
        for g = (g0(i,j)):10,
            % loop on the relationships and cell features
            for l = 1:6,
                if abs(coef(g,j,i,l))>0,
                    % concatenation of all the squared correlation
                    % coefficients
                    X = [X,(coef(g,j,i,l))^2];
                end
            end
        end
    end
end

% computing the histogram
edges = (0:0.1:1);
edges = round(edges*100)/100;
[N,BIN] = histc(X,edges);

% and plotting it
h = bar(edges,N,'histc');
delete(findobj('marker','*'));
set(gca, 'FontSize', 25, 'fontName','Times');
set(gca,'XTick',edges);
xlim([0 1])
set(h,'facecolor',[0.3 0.3 0.3]);
xlabel('R^2');
ylabel('Number of couple of distributions');

figure(12); saveas(gcf,'figures/R_squared');

end

