function [] = mainmodelstat( i,ndraw, Data, Options, pop_parameters )
%MAINMODELSTAT computes randomly generated cell lineages based on cell 
%populations statistics for embryo i - it then draws embryo level dynamics
%obtained on artificial cell lineages and compare them to the actual embryo

if i == 1,
    ng = [1,1,2,1];
    abslim = [4.25 6.3];
elseif i == 2,
    ng = [0,0,0,0];
    abslim = [6 7.3];
elseif i == 3,
    ng = [0,2,3,3];
    abslim = [4.2 9.5];
elseif i == 4,
    ng = [1,1,2,1];
    abslim = [4.85 7.8];
elseif i ==5,
    ng = [0,1,3,2];
    abslim = [4.75 7.5];
end

% generating random cell lineages
[nbcellsstatmean, volcellsstatmean, surfcellsstatmean ] = modelstat( i,ng,ndraw,Data, Options, pop_parameters);
% visualising embryo level dynamics based on these random lineages
modelstat_visualization(i,nbcellsstatmean,volcellsstatmean,surfcellsstatmean,Data, Options,abslim);

end

