function [rescaledtime,a,Nbcells] = temporalrescaling(Data, Options)
%TEMPORALRESCALING computes temporal rescaling of the embryos in the cohort
%based on the total number of cells

exp_time = zeros(max(Options.tmax),Options.Nb);
Nbcells = zeros(max(Options.tmax),5,Options.Nb);

for i = 1:Options.Nb, % for each embryo in the cohort
    % a time vector using experimental parameters
    exp_time(1:Options.tmax(i),i) = (1:Options.tmax(i))*(Options.deltat(i))+Options.tinit(i);
    
    % computing the number of cells through time for the whole embryo and
    % for each population of cells
    for k = 1:length(Data.L{i}(:,1)),
        % the total number of cells is incremented through the whole
        % cell cycle of the cell
        % for the whole embryo
        Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),5,i) = 1+Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),5,i);
        % and for the corresponding cell population
        Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),Data.L{i}(k,2),i) = 1+Nbcells(Data.L{i}(k,3):Data.L{i}(k,4),Data.L{i}(k,2),i);
    end
end

% a grid on which the function t(n) will be interpolated is constructed
i = 3;
nmini = Nbcells(Options.tmax(i),5,i);
nmaxi = Nbcells(1,5,i);
% minimal and maximal number of cells among all embryos and all time steps
for i = 1:Options.Nb,
    nmini = min(nmini, Nbcells(1,5,i));
    nmaxi = max(nmaxi, Nbcells(Options.tmax(i),5,i));
end
% the grid spans the minimal number of cells to the maximal number of cells
grid = (floor(nmini):10:ceil(nmaxi));
interpolated_exp_time = zeros(length(grid),Options.Nb);
%%
% constructing the function t(n) (time as a function of the number of
% cells, the inverse of the number cells through time N(t))
for i = 1:Options.Nb,
    xcon = [];
    ycon = [];
    ind = [];
    % we use the function consolidator to inverse the function N(t) as it
    % is not a bijection - when several time t1, t2, .. correspond to the
    % same number of cells N, the we choose the minimal one for the value
    % of the function t(N)
    [xcon,ycon,ind] = consolidator(Nbcells(1:Options.tmax(i),5,i),exp_time(1:Options.tmax(i),i),'min');
    % the function t(N) is then interpolated on the common grid for each
    % embryo
    interpolated_exp_time(:,i) = interp1(xcon,ycon,grid);
    % with 0 instead of NaN
    indtemp = isnan(interpolated_exp_time(:,i));
    interpolated_exp_time(indtemp,i) = zeros(length(find(indtemp)),1);
end

% the different t(n) of the cohort are averaged on the grid
mean_exp_time = zeros(size(grid));
for k = 1:length(grid),
    mean_exp_time(k) = mean(interpolated_exp_time(k,find(interpolated_exp_time(k,:))));
end

%%
% Optimal affine transformation parameters of the time are obtained by mean
% least square between individual t(n) and averaged t(n)

a = zeros(Options.Nb,2);
for i = 1:5,
    indtemp = find(interpolated_exp_time(:,i));
    ind = find(indtemp<= 35);
    f = @(x) sum((mean_exp_time(indtemp(ind))'-interpolated_exp_time(indtemp(ind),i)*x(1) - x(2)).^2);  % The parametrized function.
    % The parameters
    c = (fminsearch(@(x) f(x),[1,10]));
    a(i,1) = c(1);
    a(i,2) = c(2);
end

%% visualization of the number of cells function before and after rescaling
figure(1),
subplot(2,2,1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:Options.Nb,
    indtemp = find(exp_time(:,i));
    plot(exp_time(indtemp,i)/60,Nbcells(indtemp,5,i),'-','Color',Options.colors_embryos(i,:),'LineWidth',5);
    hold on,
end
xlim([220 710]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min');
ylabel('Number of cells');
title('Before temporal rescaling');

figure(1),
subplot(2,2,2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:Options.Nb,
    indtemp = find(exp_time(:,i));
    plot((a(i,1)*exp_time(indtemp,i)+a(i,2))/60,Nbcells(indtemp,5,i),'-','Color',Options.colors_embryos(i,:),'LineWidth',5);
    hold on,
end
xlim([220 710]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min');
ylabel('Number of cells');
title('After temporal rescaling');

%%
% rescaled time vector are stored
rescaledtime = cell(Options.Nb,1);

for i = 1:Options.Nb,
    indtemp = find(exp_time(:,i));
    rescaledtime{i} =  (a(i,1)*exp_time(indtemp,i)+a(i,2))';
end

%% visualization of the coefficients
figure(1),
subplot(2,2,3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i =1:Options.Nb,
    hold on
    plot(a(i,1),a(i,2)/60,'o','MarkerSize',25,'MarkerEdgeColor','none','MarkerFaceColor',Options.colors_embryos(i,:));
end
xlim([0.7 1.5]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('Scaling coefficient');
ylabel('Translation coefficient (min)');
title('Coefficients for the affine transformation');

% saving the output figure
figure(1); saveas(gcf,'figures/temporal_rescaling');

end