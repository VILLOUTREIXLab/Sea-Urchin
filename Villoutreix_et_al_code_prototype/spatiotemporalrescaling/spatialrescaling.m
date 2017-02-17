function [ Data ] = spatialrescaling(Data, Options)
%SPATIALRESCALING Computes spatial rescaling of the embryos in the cohort
%based on the total cell volume

exp_time = zeros(max(Options.tmax),Options.Nb);

%%
% spatial rescaling coefficients are initialized
a = zeros(5,1);

% total surface area and volume at each time step are computed for each embryo
Volcells = zeros(max(Options.tmax),5,Options.Nb);
Surfcells = zeros(max(Options.tmax),5,Options.Nb);

for i = 1:Options.Nb,
    % non-rescaled time vector is constructed
    exp_time(1:Options.tmax(i),i) = (1:Options.tmax(i))*(Options.deltat(i))+Options.tinit(i);
    for t = 1:Options.tmax(i),
       % total volume and surface area are added at each time step ...
       % ... for the whole embryo
       indtemp = find(Data.selection{i}(:,3) == t); 
       Volcells(t,5,i) = sum(Data.selection{i}(indtemp,6));
       Surfcells(t,5,i) = sum(Data.selection{i}(indtemp,7));
       % ... for each cell population
       for j = 1:4,
           indtemp = find(Data.selection{i}(:,3) == t & Data.selection{i}(:,2) == j);
           Volcells(t,j,i) = sum(Data.selection{i}(indtemp,6));
           Surfcells(t,j,i) = sum(Data.selection{i}(indtemp,7));
       end
    end
end

%% mean volume on a grid 
% a common time grid is constructed among the cohort of embryos
tmini = Data.rescaledtime{i}(Options.tmax(i))/60;
tmaxi = Data.rescaledtime{i}(1)/60;
for i =1:5,
   tmini = min(tmini,Data.rescaledtime{i}(1)/60); 
   tmaxi = max(tmaxi,Data.rescaledtime{i}(Options.tmax(i))/60);
end

grid = (floor(tmini):ceil(tmaxi));
% interpolated volume V(t) of each embryo on the common grid
for i = 1:5,
    xcon = [];
    ycon = [];
    ind = [];
    [xcon,ycon,ind] = consolidator((Data.rescaledtime{i}(1:Options.tmax(i))/60)',Volcells(1:Options.tmax(i),5,i));% ./volmoy108));
    Volcellsinterp(:,i) = interp1(xcon,ycon,grid);
end

Volcellsmean = zeros(length(grid),1);

% average of the individual embryo volumes V(t) on the grid
for k = 1:length(grid),
    Volcellsmean(k) = nanmean(Volcellsinterp(k,find(Volcellsinterp(k,:))),2);
end

%% optimal a parameter is obtained by least square minimization

for i = 1:5,
    indtemp = find(Volcellsinterp(:,i)>0 & Volcellsinterp(:,i) ~= NaN);
    % The parameterized function.
    f = @(x) sum((Volcellsmean(indtemp)-x*Volcellsinterp(indtemp,i)).^2);  
    % The parameter.
    a(i) = (fminsearch(@(x) f(x),[1]))^(1/3);
end



%% visualization of the volume before and after rescaling

figure(2), 
subplot(2,3,1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:5,
    indtemp = 1:Options.tmax(i);
    hold on,
    plot(exp_time(indtemp,i)/60,Volcells(indtemp,5,i),'-','Color',Options.colors_embryos(i,:),'LineWidth',5);
    hold off,
end
xlim([220 710]);
ylim([0 520000]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min pf');
ylabel('\mum^3');
title('Volume before rescaling');

figure(2), 
subplot(2,3,2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
for i = 1:5,
    indtemp = find(Volcells(:,j,i));
    plot(Data.rescaledtime{i}(indtemp)/60,a(i)^3*Volcells(indtemp,j,i),'LineWidth',5,'Color',Options.colors_embryos(i,:));
hold on,
end
end
ylim([0 520000]);
xlim([220 710]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min pf');
ylabel('\mum^3');
title('Volume after rescaling');

%% visualization of the surface area before and after rescaling
figure(2), 
subplot(2,3,4)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i = 1:5,
    indtemp = 1:Options.tmax(i);
    hold on,
    plot(exp_time(indtemp,i)/60,Surfcells(indtemp,5,i),'-','Color',Options.colors_embryos(i,:),'LineWidth',5);
    hold off,
end
xlim([220 710]);
ylim([0 280000]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min pf');
ylabel('\mum^2');
title('Surface area before rescaling');

figure(2), 
subplot(2,3,5)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
for i = 1:5,
    indtemp = find(Surfcells(:,j,i));
    plot(Data.rescaledtime{i}(indtemp)/60,a(i)^2*Surfcells(indtemp,j,i),'LineWidth',5,'Color',Options.colors_embryos(i,:));
hold on,
end
end
ylim([0 280000]);
xlim([220 710]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('min pf');
ylabel('\mum^2');
title('Surface area after rescaling');

%% visualization of rescaling coefficients

figure(2),
subplot(2,3,3),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i =1:5,
    hold on
    plot(i,a(i),'o','MarkerSize',25,'MarkerEdgeColor','none','MarkerFaceColor',Options.colors_embryos(i,:));
    hold on,
end
xlim([0.9 5.1]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('embryo');
ylabel('Scaling coefficient');

% saving the output figure
figure(2); saveas(gcf,'figures/spatial_rescaling');

%% Updating the values of cell volumes and surface area

for i = 1:Options.Nb,
    Data.L{i}(:,9) = a(i)^3.*Data.L{i}(:,9);
    Data.L{i}(:,10) = a(i)^2.*Data.L{i}(:,10);
    Data.selection{i}(:,6) = a(i)^3.*Data.selection{i}(:,6);
    Data.selection{i}(:,7) = a(i)^2.*Data.selection{i}(:,7);
end

Data.coefRedimSpatial = a;
end

