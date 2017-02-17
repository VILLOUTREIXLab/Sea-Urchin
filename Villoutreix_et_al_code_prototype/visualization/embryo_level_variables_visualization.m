function [] = embryo_level_variables_visualization(Data, Options)
%EMBRYO_LEVEL_VARIABLES_VISUALISATION provides plots of embryo level
%dynamics as a function of (rescaled) time

%% Number of cells

figure(3),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
    for i = 1:Options.Nb,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.NbCells(1:Options.tmax(i),j,i),'Color',Options.colors_embryos(i,:),'LineWidth',5);
        hold on
    end
end

abslim = [4.2 9.5];
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Cell number');

% saving the output figure
figure(3); saveas(gcf,'figures/nbcells');

%% Volume

figure(4),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
    for i = 1:Options.Nb,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.VolCells(1:Options.tmax(i),j,i),'Color',Options.colors_embryos(i,:),'LineWidth',4);
        hold on
    end
end

abslim = [4.2 9.5];
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('Time (hpf)');
xlim(abslim)
ylabel('Total cellular volume');

% saving the output figure
figure(4); saveas(gcf,'figures/cellvolume');

%% surface area

figure(5),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for j = 1:5,
    for i = 1:Options.Nb,
        plot(Data.rescaledtime{i}(1:Options.tmax(i))/3600,Data.SurfCells(1:Options.tmax(i),j,i),'Color',Options.colors_embryos(i,:),'LineWidth',4);
        hold on,
    end
end
xlim([4.2 9.5]);
ylim([0 260000]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('Time (hpf)');
ylabel('Total cellular surface area');

% saving the output figure
figure(5); saveas(gcf,'figures/cellsurface');

end

