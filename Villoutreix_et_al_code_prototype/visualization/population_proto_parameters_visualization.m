function [] = population_proto_parameters_visualization( proto_parameters, Options )
%POPULATION_PROTO_PARAMETERS_VISUALISATION Provides visualisations of population
%statistical parameters for the prototype: division time, cell cycle
%length, (log)volume,(log)surface coefficient, daughter/mother
%volume/surface area coefficient.
%The parameters are shown for each embryo, population, generation.
%
%This script is very similar to population_parameters_visualization.m -
%showing only statistics for the prototype

%%% the same layout is kept for each statistics - comments are provided for
%%% division time

%% division time
figure(16),
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
        % here only one point will be drawn in each column
        dt = 1;
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
        % if there is a prototypical statistics for the current couple
        % selection, generation
        if (proto_parameters.meandivisiontimeProto(g,j,1)>0),
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meandivisiontimeProto(g,j,1),sqrt(proto_parameters.meandivisiontimeProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
title('Division time');

%% Cell cycle length
figure(16),
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
        dt = 1;
        if j == 2,
            abscissa = 1+(g-7);
        elseif j == 1,
            abscissa = 0;
        elseif j == 3,
            abscissa = 3+(g-7);
        else
            abscissa = 6+(g-7);
        end
        if (proto_parameters.meancyclelengthProto(g,j,1)>0),
            hold on
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meancyclelengthProto(g,j,1),sqrt(proto_parameters.meancyclelengthProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
figure(16),
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
    for g = 6:10,
        c = 0;
        dt = 1;
        if j == 2,
            abscissa = 3+(g-7);
        elseif j == 1,
            abscissa = 0+(g-6);
        elseif j == 3,
            abscissa = 7.5+ (g-7-0.5);
        else
            abscissa = 12.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogvolumeProto(g,j,1)>0),
            hold on
            errorbar2(abscissa+(2*c+1)/(2*dt),exp(proto_parameters.meanlogvolumeProto(g,j,1)+0.5*proto_parameters.meanlogvolumeProto(g,j,2)),sqrt((exp(proto_parameters.meanlogvolumeProto(g,j,2))-1)*(exp(2*proto_parameters.meanlogvolumeProto(g,j,1)+proto_parameters.meanlogvolumeProto(g,j,2)))),'k');
            hold on;
            c = c+1;
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
title('Volume');


%% coefficient alpha volume
figure(18),
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
        dt = 1;
        if j == 2,
            abscissa = 1+(g-7);
        elseif j == 1,
            abscissa = 0;
        elseif j == 3,
            abscissa = 4.5+ (g-7-0.5);
        else
            abscissa = 8.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogcoefvolProto(g,j,1)~=0),
            hold on
            errorbar2(abscissa+(2*c+1)/(2*dt),exp(proto_parameters.meanlogcoefvolProto(g,j,1)+0.5*proto_parameters.meanlogcoefvolProto(g,j,2)),sqrt((exp(proto_parameters.meanlogcoefvolProto(g,j,2))-1)*(exp(2*proto_parameters.meanlogcoefvolProto(g,j,1)+proto_parameters.meanlogcoefvolProto(g,j,2)))),'k');
            hold on;
            c = c+1;
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
figure(16),
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
        dt = 1;
        if j == 2,
            abscissa = 3+(g-7);
        elseif j == 1,
            abscissa = 0+(g-6);
        elseif j == 3,
            abscissa = 7.5+ (g-7-0.5);
        else
            abscissa = 12.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogsurfaceProto(g,j,1)>0),
            hold on
            errorbar2(abscissa+(2*c+1)/(2*dt),exp(proto_parameters.meanlogsurfaceProto(g,j,1)+0.5*proto_parameters.meanlogsurfaceProto(g,j,2)),sqrt((exp(proto_parameters.meanlogsurfaceProto(g,j,2))-1)*(exp(2*proto_parameters.meanlogsurfaceProto(g,j,1)+proto_parameters.meanlogsurfaceProto(g,j,2)))),'k');
            hold on;
            c = c+1;
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
title('Statistics for the surface');

%% coefficient beta surface
figure(18),
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
    for g = 6:10,
        c = 0;
        dt = 1;
        if j == 2,
            abscissa = 1+(g-7);
        elseif j == 1,
            abscissa = 0;
        elseif j == 3,
            abscissa = 4.5+ (g-7-0.5);
        else
            abscissa = 8.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogcoefsurfProto(g,j,1)~=0),
            hold on
            errorbar2(abscissa+(2*c+1)/(2*dt),exp(proto_parameters.meanlogcoefsurfProto(g,j,1)+0.5*proto_parameters.meanlogcoefsurfProto(g,j,2)),sqrt((exp(proto_parameters.meanlogcoefsurfProto(g,j,2))-1)*(exp(2*proto_parameters.meanlogcoefsurfProto(g,j,1)+proto_parameters.meanlogcoefsurfProto(g,j,2)))),'k');
            hold on;
            c = c+1;
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
figure(17),
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
        dt = 1;
        if j == 2,
            abscissa = 3+(g-7);
        elseif j == 1,
            abscissa = 0+(g-6);
        elseif j == 3,
            abscissa = 7.5+(g-7-0.5);
        else
            abscissa = 12.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogvolumeProto(g,j,1)>0),
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meanlogvolumeProto(g,j,1),sqrt(proto_parameters.meanlogvolumeProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
ylabel('\mum^3');
title('Statistics for the log - volume');

%% log - coefficient alpha volume
figure(17)
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
        dt = 1;
        if j == 2,
            abscissa = 1+(g-7);
        elseif j == 1,
            abscissa = 0;
        elseif j == 3,
            abscissa = 4.5+(g-7-0.5);
        else
            abscissa = 8.5+(g-7-0.5);
        end
        if (abs(proto_parameters.meanlogcoefvolProto(g,j,1))>0),
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meanlogcoefvolProto(g,j,1),sqrt(proto_parameters.meanlogcoefvolProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
figure(17),
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
        dt = 1;
        if j == 2,
            abscissa = 3+(g-7);
        elseif j == 1,
            abscissa = 0+(g-6);
        elseif j == 3,
            abscissa = 7.5+(g-7-0.5);
        else
            abscissa = 12.5+(g-7-0.5);
        end
        if (proto_parameters.meanlogsurfaceProto(g,j,1)>0),
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meanlogsurfaceProto(g,j,1),sqrt(proto_parameters.meanlogsurfaceProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
ylabel('\mum^2');
title('Statistics for the log - surface');

%% LOG - coefficient beta surface
figure(17)
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
        dt = 1;
        if j == 2,
            abscissa = 1+(g-7);
        elseif j == 1,
            abscissa = 0;
        elseif j == 3,
            abscissa = 4.5+(g-7-0.5);
        else
            abscissa = 8.5+(g-7-0.5);
        end
        if (abs(proto_parameters.meanlogcoefsurfProto(g,j,1))>0),
            errorbar2(abscissa+(2*c+1)/(2*dt),proto_parameters.meanlogcoefsurfProto(g,j,1),sqrt(proto_parameters.meanlogcoefsurfProto(g,j,2)),'k');
            hold on;
            c = c+1;
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
ylabel('Min pf');
title('Statistics for log B coefficient');

%% saving figures
figure(16); saveas(gcf,'figures/proto_parameters_1');
figure(17); saveas(gcf,'figures/proto_parameters_2');
figure(18); saveas(gcf,'figures/proto_parameters_3');

end

