function [] = population_parameters_visualization( pop_parameters, Options )
%POPULATION_PARAMETERS_VISUALISATION Provides visualisations of population
%statistical parameters: division time, cell cycle length, (log)volume,
%(log)surface coefficient, daughter/mother volume/surface area coefficient.
%The parameters are shown for each embryo, population, generation

I = 1:Options.Nb; %cohort of embryos

%%% the same layout is kept for each statistics - comments are provided for
%%% division time

%% Division time
figure(6),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,2,1)

% the layout is divided in columns corresponding to populations and
% generations

% the first loop is on populations
for j = 1:4,
    % the background is filled with a color corresponding to the cell
    % population
    if j == 1,
        hold on,
        fill([0,1,1,0],[200,200,600,600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,4,4,1],[200,200,600,600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([4,8,8,4],[200,200,600,600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([8,12,12,8],[200,200,600,600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    
    % the second loop is on the generations
    for g = 6:9
        c = 0;
        % the number of embryos providing statistics for the given cell
        % population (defined with population and generation) is computed
        % as dt and will define the horizontal distribution in the
        % considered column
        dt = length(find(pop_parameters.meandivisiontime(g,j,1,:)));
        
        % the third loop is on the embryos and will define where to draw
        % points in the column population x generation
        for i = I,
            % the left abscissa of the column is defined according to the
            % couple population, generation
            if j == 2,
                abscissa = 2+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 5.5+ (g-7-0.5);
            else
                abscissa = 9.5+(g-7-0.5);
            end
            % if the set of points to be drawn is not empty
            if (pop_parameters.meandivisiontime(g,j,1,i)>0), 
                % the points are drawn sequentially among the cohort
                % the points are mean value with errorbar corresponding to
                % standard deviation
                % points are colored according to embryo
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meandivisiontime(g,j,1,i),sqrt(pop_parameters.meandivisiontime(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        % a separator is drawn between the columns
        plot(abscissa*ones(401,1),200:600,'k--');
        alpha(0.1);
    end
end
% dimensions of the figure are fixed
xlim([0 12])
ylim([260 570])
% font
set(gca, 'FontSize', 20, 'fontName','Times');
% choice of label for the axes
set(gca,'Box','off');
set(gca,'XTick',[0.5 2.5 6 10])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('Min pf');
title('Moment of division');

%% cell cycle length
figure(6),
subplot(2,2,2)
for j = 1:4,
    if j == 1,
        hold on,
        fill([0,1,1,0],[40,40,180,180],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,3,3,1],[40,40,180,180],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([3,6,6,3],[40,40,180,180],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([6,9,9,6],[40,40,180,180],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 7:9,
        c = 0;
        dt = length(find(pop_parameters.meancyclelength(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 1+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 3+(g-7);
            else
                abscissa = 6+(g-7);
            end
            if (pop_parameters.meancyclelength(g,j,1,i)>0),
                hold on
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meancyclelength(g,j,1,i),sqrt(pop_parameters.meancyclelength(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end 
        end
        plot(abscissa*ones(141,1),40:180,'k--');
        alpha(0.1);
    end
end
xlim([0 9]);
ylim([40 180]);
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[0.5 2 4.5 7.5])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('Min');
title('Cell cycle length');


%% volume
figure(6),
subplot(2,2,3)
for j = 1:4,
    if j == 1,
        hold on,
        fill([0,2,2,0],[250,250,18600,18600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([2,6,6,2],[250,250,18600,18600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([6,11,11,6],[250,250,18600,18600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([11,16,16,11],[250,250,18600,18600],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10
        c = 0;
        dt = length(find(pop_parameters.meanvolume(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 3+(g-7);
            elseif j == 1,
                abscissa = 0+(g-6);
            elseif j == 3,
                abscissa = 7.5+ (g-7-0.5);
            else
                abscissa = 12.5+(g-7-0.5);
            end
            if (pop_parameters.meanvolume(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meanvolume(g,j,1,i),sqrt(pop_parameters.meanvolume(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        plot(abscissa*ones(18351,1),250:18600,'k--');
        alpha(0.1);
    end
end
xlim([0 16])
ylim([250 18600])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[1 4 8.5 13.5])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('\mum^3');
title('Statistics for the volume');


%% coefficient alpha volume
figure(8),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,1,1),

for j = 1:4,
    if j == 1,
        hold on,
        fill([0,1,1,0],[0.2,0.2,0.85,0.85],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,4,4,1],[0.2,0.2,0.85,0.85],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([4,8,8,4],[0.2,0.2,0.85,0.85],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([8,12,12,8],[0.2,0.2,0.85,0.85],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10
        c = 0;
        dt = length(find(pop_parameters.meancoefvol(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 1+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 4.5+ (g-7-0.5);
            else
                abscissa = 8.5+(g-7-0.5);
            end
            if (pop_parameters.meancoefvol(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meancoefvol(g,j,1,i),sqrt(pop_parameters.meancoefvol(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;      
            end
        end
        plot(abscissa*ones(66,1),0.2:0.01:0.85,'k--');
    end
end
xlim([0 12])
ylim([0.2 0.85])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[0.5 2.5 6 10])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
title('Statistics for the A coefficient');


%% surface
figure(6),
subplot(2,2,4),

for j = 1:4,
    if j == 1,
        hold on,
        fill([0,2,2,0],[200,200,4000,4000],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([2,6,6,2],[200,200,4000,4000],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([6,11,11,6],[200,200,4000,4000],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([11,16,16,11],[200,200,4000,4000],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10,
        c = 0;
        dt = length(find(pop_parameters.meansurface(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 3+(g-7);
            elseif j == 1,
                abscissa = 0+(g-6);
            elseif j == 3,
                abscissa = 7.5+ (g-7-0.5);
            else
                abscissa = 12.5+(g-7-0.5);
            end
            if (pop_parameters.meansurface(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meansurface(g,j,1,i),sqrt(pop_parameters.meansurface(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end  
        end
        plot(abscissa*ones(3801,1),200:4000,'k--');
        alpha(0.1);
    end
end
xlim([0 16])
ylim([200 4000])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[1 4 8.5 13.5])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('\mum^2');
title('Statistics for the surface area');

%% coefficient beta surface
figure(8),
subplot(2,1,2),
for j = 1:4,
    if j == 1,
        hold on,
        fill([0,1,1,0],[0.3,0.3,0.92,0.92],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,4,4,1],[0.3,0.3,0.92,0.92],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([4,8,8,4],[0.3,0.3,0.92,0.92],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([8,12,12,8],[0.3,0.3,0.92,0.92],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10
        c = 0;
        dt = length(find(pop_parameters.meancoefvol(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 1+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 4.5+ (g-7-0.5);
            else
                abscissa = 8.5+(g-7-0.5);
            end
            if (pop_parameters.meancoefsurf(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meancoefsurf(g,j,1,i),sqrt(pop_parameters.meancoefsurf(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        plot(abscissa*ones(63,1),0.3:0.01:0.92,'k--');
    end
end
xlim([0 12])
ylim([0.3 0.92])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[0.5 2.5 6 10])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
title('Statistics for the B coefficient');


%% LOG volume
figure(7),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,2,1),

for j = 1:4,
    if j == 1,
        hold on,
        fill([0,2,2,0],[6,6,10,10],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([2,6,6,2],[6,6,10,10],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([6,11,11,6],[6,6,10,10],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([11,16,16,11],[6,6,10,10],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10
        c = 0;
        dt = length(find(pop_parameters.meanlogvolume(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 3+(g-7);
            elseif j == 1,
                abscissa = 0+(g-6);
            elseif j == 3,
                abscissa = 7.5+ (g-7-0.5);
            else
                abscissa = 12.5+(g-7-0.5);
            end
            if (pop_parameters.meanlogvolume(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meanlogvolume(g,j,1,i),sqrt(pop_parameters.meanlogvolume(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        plot(abscissa*ones(length(6:0.1:10),1), 6:0.1:10,'k--');
        alpha(0.1);
    end
end
xlim([0 16])
ylim([6 10])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[1 4 8.5 13.5])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('log(\mum^3)');
title('Statistics for the log - volume');


%% log - coefficient alpha volume
figure(7)
subplot(2,2,3),

for j = 1:4,
    if j == 1,
        hold on,
        fill([0,1,1,0],[-1.4,-1.4,-0.3,-0.3],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,4,4,1],[-1.4,-1.4,-0.3,-0.3],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([4,8,8,4],[-1.4,-1.4,-0.3,-0.3],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([8,12,12,8],[-1.4,-1.4,-0.3,-0.3],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10
        c = 0;
        dt = length(find(pop_parameters.meanlogcoefvol(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 1+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 4.5+ (g-7-0.5);
            else
                abscissa = 8.5+(g-7-0.5);
            end
            if (pop_parameters.meanlogcoefvol(g,j,1,i)~=0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meanlogcoefvol(g,j,1,i),sqrt(pop_parameters.meanlogcoefvol(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        plot(abscissa*ones(length(-1.4:0.1:-0.3),1),-1.4:0.1:-0.3,'k--');
        alpha(0.1);
    end
end
xlim([0 12])
ylim([-1.4 -0.3])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[0.5 2.5 6 10])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
title('Statistics for log A coefficient');


%% LOG - surface
figure(7),
subplot(2,2,2),
for j = 1:4,
    if j == 1,
        hold on,
        fill([0,2,2,0],[5.5,5.5,8.5,8.5],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([2,6,6,2],[5.5,5.5,8.5,8.5],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([6,11,11,6],[5.5,5.5,8.5,8.5],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([11,16,16,11],[5.5,5.5,8.5,8.5],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10,
        c = 0;
        dt = length(find(pop_parameters.meanlogsurface(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 3+(g-7);
            elseif j == 1,
                abscissa = 0+(g-6);
            elseif j == 3,
                abscissa = 7.5+ (g-7-0.5);
            else
                abscissa = 12.5+(g-7-0.5);
            end
            if (pop_parameters.meanlogsurface(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meanlogsurface(g,j,1,i),sqrt(pop_parameters.meanlogsurface(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end 
        end
        plot(abscissa*ones(length(5.5:0.1:8.5),1),5.5:0.1:8.5,'k--');
        alpha(0.1);
    end
end
xlim([0 16])
ylim([5.5 8.5])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[1 4 8.5 13.5])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
ylabel('log(\mum^2)');
title('Statistics for the log - surface area');

%% LOG - coefficient beta surface
figure(7)
subplot(2,2,4),

for j = 1:4,
    if j == 1,
        hold on,
        fill([0,1,1,0],[-1,-1,-0.2,-0.2],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 2,
        hold on,
        fill([1,4,4,1],[-1,-1,-0.2,-0.2],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    elseif j == 3,
        hold on,
        fill([4,8,8,4],[-1,-1,-0.2,-0.2],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    else
        fill([8,12,12,8],[-1,-1,-0.2,-0.2],Options.colors_selections(j,:),'EdgeColor','none');
        alpha(0.1);
    end
    for g = 6:10,
        c = 0;
        dt = length(find(pop_parameters.meanlogcoefsurf(g,j,1,:)));
        for i = I,
            if j == 2,
                abscissa = 1+(g-7);
            elseif j == 1,
                abscissa = 0;
            elseif j == 3,
                abscissa = 4.5+ (g-7-0.5);
            else
                abscissa = 8.5+(g-7-0.5);
            end
            if (pop_parameters.meancoefsurf(g,j,1,i)>0),
                errorbar2(abscissa+(2*c+1)/(2*dt),pop_parameters.meanlogcoefsurf(g,j,1,i),sqrt(pop_parameters.meanlogcoefsurf(g,j,2,i)),Options.colors_embryos(i,:));
                hold on;
                c = c+1;
            end
        end
        plot(abscissa*ones(length(-1:0.1:0.2),1),-1:0.1:0.2,'k--');
        alpha(0.1);
    end
end
xlim([0 12])
ylim([-1 -0.2])
set(gca, 'FontSize', 20, 'fontName','Times');
set(gca,'Box','off');
set(gca,'XTick',[0.5 2.5 6 10])
set(gca,'XTickLabel',{'SMic';'LMic';'Mac';'Mes'});
title('Statistics for log B coefficient');



%%

% saving the output figure
figure(6); saveas(gcf,'figures/cell_population_parameters_1');
figure(7); saveas(gcf,'figures/cell_population_parameters_2');
figure(8); saveas(gcf,'figures/cell_population_parameters_3');

end

