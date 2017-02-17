function [] = modelevaluation( pop_parameters,Options )
%MODELEVALUATION systematically compares the analytical results of the
%model with the empirical results

%%% model evaluation with Kullback Leibler divergence
% goal: compare mean and standard deviation of empirical statistics with
% those obtained from the model as a sum of independent variables

% initialization of variables
% symmetrized-kullback-leibler divergence
val = zeros(4,5,3);
% score defined as a ratio val/val2
val1 = zeros(4,5,3);
% mean distance between embryo and the cohort
val2 = zeros(4,5,3);

% cohort
I = 1:Options.Nb;


%---------------
% division time
%---------------

% loop over each embryo
for i = I,
    % loop over the populations
    for j = 1:4,
        %  find the values of mean division time at all generations
        indtemp = find(pop_parameters.meandivisiontime(:,j,1,i));
        if ~isempty(indtemp),
            % find the generation corresponding to the last cell cycle
            gfin = indtemp(end);
            if length(find(pop_parameters.meandivisiontime(gfin,j,1,:))) == 1,
                gfin = gfin - 1;
            end
            % find the generation corresponding to the initial cell cycle
            ginit = indtemp(1);
            if gfin > ginit,
                % empirical values of the mean and standard deviation of
                % division time for the last cell cycle
                meanemp = pop_parameters.meandivisiontime(gfin,j,1,i);
                varemp = pop_parameters.meandivisiontime(gfin,j,2,i);
                
                % the mean value obtained from the model is the sum of the
                % previous generations mean cell cycle lengths and initial
                % mean division time
                meanmodele = sum(pop_parameters.meancyclelength((ginit+1:gfin),j,1,i)) + pop_parameters.meandivisiontime(ginit,j,1,i);
                % similarly for the variance
                varmodele = sum(pop_parameters.meancyclelength((ginit+1:gfin),j,2,i)) + pop_parameters.meandivisiontime(ginit,j,2,i);
                
                val(j,i,1) = 0;
                % we can now compute the Kullback-Leibler symmetrized
                % divergence assuming that the distributions are normally
                % distributed
                val(j,i,1) = 0.25*(varmodele/varemp + varemp/varmodele + (meanmodele - meanemp)^2*(1/varemp + 1/varmodele) - 2);
                
                % compare to other embryos so that the difference between
                % the model and the data can be normalized with respect to
                % the spread over the cohort
                
                % find the values of the mean division time among all the
                % embryos of the cohort for the last generation
                indtemp = find(pop_parameters.meandivisiontime(gfin,j,1,:));
                
                if length(indtemp)>1,
                    c = 0;
                    val2(j,i,1) = 0;
                    % for each embryo of the cohort
                    for k = 1:length(indtemp),
                        if indtemp(k) ~= i,
                            % compute the symmetrized kullback-leibler
                            % divergence between embryo k and embryo i
                            restemp = 0.25*(pop_parameters.meandivisiontime(gfin,j,2,indtemp(k))/pop_parameters.meandivisiontime(gfin,j,2,i) + pop_parameters.meandivisiontime(gfin,j,2,i)/pop_parameters.meandivisiontime(gfin,j,2,indtemp(k)) + (pop_parameters.meandivisiontime(gfin,j,1,i) - pop_parameters.meandivisiontime(gfin,j,1,indtemp(k)))^2*(1/pop_parameters.meandivisiontime(gfin,j,2,indtemp(k)) + 1/pop_parameters.meandivisiontime(gfin,j,2,i))-2);
                            val2(j,i,1) = val2(j,i,1)+restemp;
                            c = c+1;
                        end
                    end
                    % average the difference with the other embryos
                    val2(j,i,1) = val2(j,i,1)/c;
                end
                if val(j,i,1)>0,
                    if val2(j,i,1)>0,
                        % finally, the result that will be used is the
                        % difference between the measure and the model
                        % normalized by the difference to other embryos
                        val1(j,i,1) = val(j,i,1)/val2(j,i,1);
                    end
                end
            end
        end
    end
end

% putting all model evaluations concerning division time together
res_divisiontime = [];
for i = I,
    for j = 1:4,
        if val1(j,i,1)>0,
            res_divisiontime = [res_divisiontime,val1(j,i,1)];
        end
    end
end

%---------------
% volume
%---------------

% loop over each embryo
for i = I,
    % loop over the populations
    for j = 1:4,
        %  find the values of mean log(volume) at all generations
        indtemp = find(pop_parameters.meanlogvolume(:,j,1,i));
        if ~isempty(indtemp),
            % find the generation corresponding to the last cell cycle
            gfin = indtemp(end);
            if length(find(pop_parameters.meanlogvolume(gfin,j,1,:))) == 1,
                gfin = gfin - 1;
            end
            % find the generation corresponding to the initial cell cycle
            ginit = indtemp(1);
            if gfin > ginit,
                % empirical values of the mean and standard deviation of
                % log volume for the last cell cycle
                meanemp = pop_parameters.meanlogvolume(gfin,j,1,i);
                varemp = pop_parameters.meanlogvolume(gfin,j,2,i);
                
                % the mean value obtained from the model is the sum of the
                % previous generations mean log daughter/mother volume
                % ratio and initial mean log volume
                meanmodele = sum(pop_parameters.meanlogcoefvol((ginit+1:gfin),j,1,i)) + pop_parameters.meanlogvolume(ginit,j,1,i);
                % similarly for the variance
                varmodele = sum(pop_parameters.meanlogcoefvol((ginit+1:gfin),j,2,i)) + pop_parameters.meanlogvolume(ginit,j,2,i);
                
                val(j,i,2) = 0;
                % we can now compute the Kullback-Leibler symmetrized
                % divergence assuming that the distributions are normally
                % distributed
                val(j,i,2) = 0.25*(varmodele/varemp + varemp/varmodele + (meanmodele - meanemp)^2*(1/varemp + 1/varmodele) - 2);
                
                % compare to other embryos so that the difference between
                % the model and the data can be normalized with respect to
                % the spread over the cohort
                
                % find the values of the mean log volume among all the
                % embryos of the cohort for the last generation
                indtemp = find(pop_parameters.meanlogvolume(gfin,j,1,:));
                
                if length(indtemp)>1,
                    val2(j,i,2) = 0;
                    c = 0;
                    % for each embryo of the cohort
                    for k = 1:length(indtemp),
                        if indtemp(k) ~= i,
                            % compute the symmetrized kullback-leibler
                            % divergence between embryo k and embryo i
                            restemp = 0.25*(pop_parameters.meanlogvolume(gfin,j,2,indtemp(k))/pop_parameters.meanlogvolume(gfin,j,2,i) + pop_parameters.meanlogvolume(gfin,j,2,i)/pop_parameters.meanlogvolume(gfin,j,2,indtemp(k)) + (pop_parameters.meanlogvolume(gfin,j,1,i) - pop_parameters.meanlogvolume(gfin,j,1,indtemp(k)))^2*(1/pop_parameters.meanlogvolume(gfin,j,2,indtemp(k)) + 1/pop_parameters.meanlogvolume(gfin,j,2,i))-2);
                            val2(j,i,2) = val2(j,i,2)+restemp;
                            c = c+1;
                        end
                    end
                    % average the difference with the other embryos
                    val2(j,i,2) = val2(j,i,2)/c;
                end
                if val(j,i,2)>0,
                    if val2(j,i,2)>0,
                        % finally, the result that will be used is the
                        % difference between the measure and the model
                        % normalized by the difference to other embryos
                        val1(j,i,2) = val(j,i,2)/val2(j,i,2);
                    end
                end
            end
        end
    end
end

% putting all model evaluations concerning volume
res_volume = [];

for i = I,
    for j = 1:4,
        if val1(j,i,2)>0,
            res_volume = [res_volume,val1(j,i,2)];
        end
    end
end

%---------------
% surface area
%---------------

% loop over each embryo
for i = I,
    % loop over the populations
    for j = 1:4,
        %  find the values of mean log(surface area) at all generations
        indtemp = find(pop_parameters.meanlogsurface(:,j,1,i));
        if ~isempty(indtemp),
            % find the generation corresponding to the last cell cycle
            gfin = indtemp(end);
            if length(find(pop_parameters.meanlogsurface(gfin,j,1,:))) == 1,
                gfin = gfin - 1;
            end
            % find the generation corresponding to the initial cell cycle
            ginit = indtemp(1);
            if gfin > ginit,
                % empirical values of the mean and standard deviation of
                % log surface area for the last cell cycle
                meanemp = pop_parameters.meanlogsurface(gfin,j,1,i);
                varemp = pop_parameters.meanlogsurface(gfin,j,2,i);
                
                % the mean value obtained from the model is the sum of the
                % previous generations mean log daughter/mother surface area
                % ratio and initial mean log surface area
                meanmodele = sum(pop_parameters.meanlogcoefsurf((ginit+1:gfin),j,1,i)) + pop_parameters.meanlogsurface(ginit,j,1,i);
                % similarly for the variance
                varmodele = sum(pop_parameters.meanlogcoefsurf((ginit+1:gfin),j,2,i)) + pop_parameters.meanlogsurface(ginit,j,2,i);
                val(j,i,3) = 0;
                
                % we can now compute the Kullback-Leibler symmetrized
                % divergence assuming that the distributions are normally
                % distributed
                val(j,i,3) = 0.25*(varmodele/varemp + varemp/varmodele + (meanmodele - meanemp)^2*(1/varemp + 1/varmodele) - 2);
                
                % compare to other embryos so that the difference between
                % the model and the data can be normalized with respect to
                % the spread over the cohort
                
                % find the values of the mean log surface area among all
                % the embryos of the cohort for the last generation
                indtemp = find(pop_parameters.meanlogsurface(gfin,j,1,:));
                
                if length(indtemp)>1,
                    val2(j,i,3) = 0;
                    c = 0;
                    % for each embryo of the cohort
                    for k = 1:length(indtemp),
                        if indtemp(k) ~= i,
                            % compute the symmetrized kullback-leibler
                            % divergence between embryo k and embryo i
                            restemp = 0.25*(pop_parameters.meanlogsurface(gfin,j,2,indtemp(k))/pop_parameters.meanlogsurface(gfin,j,2,i) + pop_parameters.meanlogsurface(gfin,j,2,i)/pop_parameters.meanlogsurface(gfin,j,2,indtemp(k)) + (pop_parameters.meanlogsurface(gfin,j,1,i) - pop_parameters.meanlogsurface(gfin,j,1,indtemp(k)))^2*(1/pop_parameters.meanlogsurface(gfin,j,2,indtemp(k)) + 1/pop_parameters.meanlogsurface(gfin,j,2,i))-2);
                            val2(j,i,3) = val2(j,i,3)+restemp;
                            c = c+1;
                        end
                    end
                    % average the difference with the other embryos
                    val2(j,i,3) = val2(j,i,3)/c;
                end
                if val(j,i,3)>0,
                    if val2(j,i,3)>0,
                        % finally, the result that will be used is the
                        % difference between the measure and the model
                        % normalized by the difference to other embryos
                        val1(j,i,3) = val(j,i,3)/val2(j,i,3);
                    end
                end
            end
        end
    end
end

% putting all model evaluations concerning surface area together
res_surfacearea = [];
c = 1;
for i = I,
    for j = 1:4,
        if val1(j,i,3)>0,
            res_surfacearea = [res_surfacearea,val1(j,i,3)];
            c = c+1;
        end
    end
end

%---------------
% Results visualization
%---------------

figure(15),
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% all the evaluation scores are stored in a common vector X
X = [res_divisiontime,res_volume,res_surfacearea];

% compute and plot the histogram of X
edges = (0:0.5:12);
[N,BIN] = histc(X,edges);
h = bar(edges,N,'histc');
delete(findobj('marker','*'));
set(gca,'XTick',edges);

xlim([0 11.5]);
ylim([0 37]);

set(h,'facecolor',[0.3 0.3 0.3]);
set(gca, 'FontSize', 25, 'fontName','Times');
xlabel('Normalized symmetrized Kullback Leibler Divergence between the data and the model');
ylabel('Number of distributions');
title('All statistics');

figure(15); saveas(gcf,'figures/model_evaluation');
close(figure(15));

end

