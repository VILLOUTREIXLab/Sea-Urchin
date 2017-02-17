function [] = microdynamics( Data, Options )
%MICRODYNAMICS computes averaged microdynamics over groups of cells within
%each embryo and an average over the whole cohort

% cell cycles are rescaled on a time scale of 100
nts = 100;

% volume
figure(13),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% loop on the populations (LMic, Mac, Mes)
for sp = 2:4,
    % loop on the generations
    for g = 7:9,
        % we want to draw one graph for each subpopulation and each
        % generation
        subplot(3,3,(4-sp)*3+1+g-7)
        % observation of LMic is shifted to generation 6, 7 and 8
        if sp == 2,
            g = g-1;
        end
        
        X = [];
        for i = 1:Options.Nb,
            % Computing averaged micro volume dynamics
            [m,~] = volumemicrodynamics(sp,g,i,nts,Data.L,Data.selection, Options);
            X = [X,m];
            % and plotting it
            if ~isempty(m)
                plot(1:nts,m,'Color', Options.colors_embryos(i,:));
            end
            hold on,
        end
        if ~isempty(X),
            % ploting the average over the cohort
            plot(1:nts, mean(X,2),'LineWidth',2,'Color', Options.colors_generations(g-5,:));
        end
        ylim([0.8 1.3]);
        set(gca, 'FontSize', 20, 'fontName','Times');
        title([Options.ident{sp},' gen ', num2str(g)])
        ylabel('Normalized volume')
        xlabel('Normalized cycle length')
    end
end

figure(13); saveas(gcf,'figures/volume_micro_dynamics');

% surface
figure(14),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% loop on the populations (LMic, Mac, Mes)
for sp = 2:4,
    % loop on the generations
    for g = 7:9,
        % we want to draw one graph for each subpopulation and each
        % generation
        subplot(3,3,(4-sp)*3+1+g-7)
        % observation of LMic is shifted to generation 6, 7 and 8
        if sp == 2,
            g = g-1;
        end
        X = [];
        for i = 1:Options.Nb
            % Computing averaged micro surface dynamics
            [m,] = surfacemicrodynamics(sp,g,i,nts,Data.L,Data.selection, Options);
            X = [X,m];
            % and plotting it
            if ~isempty(m)
                plot(1:nts,m,'Color', Options.colors_embryos(i,:));
            end
            hold on,
        end
        if ~isempty(X),
            % ploting the average over the cohort
            plot(1:nts, mean(X,2),'LineWidth',2,'Color', Options.colors_generations(g-5,:))
        end
        ylim([0.8 1.2]);
        set(gca, 'FontSize', 20, 'fontName','Times');
        title([Options.ident{sp},' gen ', num2str(g)])
        ylabel('Normalized surface area')
        xlabel('Normalized cycle length')
    end
end

figure(14); saveas(gcf,'figures/surface_micro_dynamics');

end

