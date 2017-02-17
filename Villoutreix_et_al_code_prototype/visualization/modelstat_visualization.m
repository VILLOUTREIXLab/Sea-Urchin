function [] = modelstat_visualization( i,nbcellsstatmean,volcellsstatmean,surfcellsstatmean,Data, Options,abslim)
%MODELSTAT_VISUALIZATION generates visualizations of embryo-level dynamics
%obtained from ensembles of artificial cell lineages and compared to
%original data

% initialization of the time grid
tinit = Data.rescaledtime{i}(1)/60;
tfin = Data.rescaledtime{i}(Options.tmax(i))/60;
time = (tinit:tfin);
ech = 1:length(time);
time = time/60;

% domain used for drawing embryo-level features
X = [time(ech)';flipud(time(ech)')];

% the embryo-level dynamics are represented with the same layout - we
% provide comments for the first one, the two others are completely similar

% ------------------------
% Number of cells in time
% ------------------------

figure(19), 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

X = [time(ech)';flipud(time(ech)')];
% for each generation and for the whole embryo
for j = 1:5,
    hold on,
    % we first draw the results of the model
    % we draw a gray area corresponding to the standard deviation of the
    % number of cells 
    % above the averaged total number of cells ...
    Y = [nbcellsstatmean(ech,j,1);flipud(nbcellsstatmean(ech,j,1)+nbcellsstatmean(ech,j,2))];
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    % ... and below the averaged total number of cells
    Y = [nbcellsstatmean(ech,j,1);flipud(nbcellsstatmean(ech,j,1)-nbcellsstatmean(ech,j,2))];    
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    
    % and we plot the averaged total number of cells
    plot(time(ech),nbcellsstatmean(ech,j,1),'-','LineWidth',5,'Color',[0.4,0.4,0.4]);
end
% then, the number of cells measured in the actual embryo in the whole 
% embryo and in each population is drawn
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),5,i),'-','LineWidth',4,'Color',Options.colors_embryos(i,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),4,i),'-','LineWidth',4,'Color',Options.colors_selections(4,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),3,i),'-','LineWidth',4,'Color',Options.colors_selections(3,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),2,i),'-','LineWidth',4,'Color',Options.colors_selections(2,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),1,i),'-','LineWidth',4,'Color',Options.colors_selections(1,:));

set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Cell number');
title('Cells number dynamics for the whole embryo and for each morphogenetic field');

figure(19); saveas(gcf,'figures/nbcells_model');

% ------------------------
% Volume in time
% ------------------------

figure(20),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
    hold on,
    X = [time(ech)';flipud(time(ech)')];
    Y = [volcellsstatmean(ech,j,1);flipud(volcellsstatmean(ech,j,1)+volcellsstatmean(ech,j,2))];
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    X = [time(ech)';flipud(time(ech)')];
    Y = [volcellsstatmean(ech,j,1);flipud(volcellsstatmean(ech,j,1)-volcellsstatmean(ech,j,2))];
    %
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    %
    
    plot(time(ech),volcellsstatmean(ech,j,1),'-','LineWidth',5,'Color',[0.4,0.4,0.4]);
end

plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),5,i),'-','LineWidth',4,'Color',Options.colors_embryos(i,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),4,i),'-','LineWidth',4,'Color',Options.colors_selections(4,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),3,i),'-','LineWidth',4,'Color',Options.colors_selections(3,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),2,i),'-','LineWidth',4,'Color',Options.colors_selections(2,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),1,i),'-','LineWidth',4,'Color',Options.colors_selections(1,:));

set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Volume in \mum^3');
title('Cells volume dynamics for the whole embryo and for each morphogenetic field');

figure(20); saveas(gcf,'figures/cellvolume_model');

% ------------------------
% Surface area in time
% ------------------------

figure(21),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for j = 1:5,
    hold on,
    X = [time(ech)';flipud(time(ech)')];
    Y = [surfcellsstatmean(ech,j,1);flipud(surfcellsstatmean(ech,j,1)+surfcellsstatmean(ech,j,2))];
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    X = [time(ech)';flipud(time(ech)')];
    Y = [surfcellsstatmean(ech,j,1);flipud(surfcellsstatmean(ech,j,1)-surfcellsstatmean(ech,j,2))];
    %
    fill(X,Y,[0.8 0.8 0.8],'EdgeColor','none');
    %
    
    plot(time(ech),surfcellsstatmean(ech,j,1),'-','LineWidth',5,'Color',[0.4,0.4,0.4]);
end

plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),5,i),'-','LineWidth',4,'Color',Options.colors_embryos(i,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),4,i),'-','LineWidth',4,'Color',Options.colors_selections(4,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),3,i),'-','LineWidth',4,'Color',Options.colors_selections(3,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),2,i),'-','LineWidth',4,'Color',Options.colors_selections(2,:));
hold on,
plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),1,i),'-','LineWidth',4,'Color',Options.colors_selections(1,:));

set(gca, 'FontSize', 20, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Surface in \mum^2');
title('Cells surface dynamics for the whole embryo and for each morphogenetic field');

figure(21); saveas(gcf,'figures/cellsurface_model');

end

