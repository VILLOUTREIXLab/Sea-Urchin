%% add paths

addpath('artificiallineages')
addpath('spatiotemporalrescaling')
addpath('statisticalanalysis')
addpath('visualization')
addpath('tools')
if ~exist('figures','dir'), 
    mkdir('figures'),
end

%% importing experimental data
clear, close all,

% the information concerning cell lineage and cell features are stored in
% two types of variables ..

% .. a selection, containing information about cells at each time step
% the fields are the following: 1) id cell 2) cell population 3) time step
% 4) id predecessor 5) global id 6) volume 7) surface area
Data.selection = {};

% .. a list of cells, containing information about cells considered
% throughout their entire cell cycle: 1) global id cell 2) cell population
% 3) initial time step 4) final time step 5) observed initial mitosis 6)
% observed final mitosis 7) generation 8) id mother cell 9) mean volume 10)
% mean surface area
Data.L = {};

for i = 1:5,
    Data.selection{i} = csvread(['dataset/selection',num2str(i),'.csv']);
    Data.L{i} = csvread(['dataset/L',num2str(i),'.csv']);
end

% Defining experimental parameter values
% Size of the cohort
Options.Nb = 5;
% Maximal time step considered for each embryo
Options.tmax = [81; 86; 90; 75; 126];
% Experimental time step in seconds
Options.deltat = [119 ; 133 ; 207 ; 220 ; 180 ];
% Experimental initial time of each experiment
Options.tinit = [13500 ; 20700 ; 15600; 18000; 15600 ];

% color code
Options.colors_selections = [[92 0 77];[255 121 255];[228 19 0];[13 173 209]]./255;
Options.colors_embryos = [[35,91,190];[39,173,25];[238,110,33];[76,77,79];[140,129,100]]./255;
Options.colors_generations = [[37,253,233];[119,181,254];[0,0,255];[102,0,153];[253,108,158];[205, 205, 13]]./255;

% names of the various subpopulations
Options.ident = {'Smic';'LMic';'Mac';'Mes'};

% number of cell in each subpopulation at generation 6
Options.ngen6 = [4,4,8,16];

% initial generation number in each population of each embryo
Options.g0 = [[6,6,6,6];[6,7,8,8];[6,6,6,6];[6,6,7,7];[6,6,7,7];];

%% Computing temporal rescaling
[Data.rescaledtime,~,~] = temporalrescaling(Data,Options);

close(figure(1));

%% Computing spatial rescaling
[Data] = spatialrescaling(Data,Options);

close(figure(2));

%% Computing and visualizing embryo level variables
[Data.NbCells, Data.VolCells, Data.SurfCells] = embryo_level_variables(Data, Options);
embryo_level_variables_visualization(Data, Options);

close(figure(3), figure(4), figure(5));

%% Computing and visualizing population level statistical parameters
[ pop_parameters ] = population_parameters( Data, Options );
population_parameters_visualization(pop_parameters, Options);

close(figure(6), figure(7), figure(8));

%% Statistical characteristics of cell feature distributions
% compute and visualize p-values associated with gaussian distribution
% goodness of fit
gaussianFit(Data, Options);
clc;
close(figure(9), figure(10), figure(11));

% Compute and visualize R-square for correlations between sisters
% and mother daughter relationships
correlations(Data, Options);

close(figure(12));

%% Micro dynamics
microdynamics(Data, Options);

close(figure(13),figure(14));

%% Model evaluation
modelevaluation(pop_parameters,Options);

close(figure(15));

%% Computing and visualizing Prototypical Population level parameters
% computing and visualizing prototypical coefficients for cell groups
% statistics
[proto_parameters] = prototypical_parameters(pop_parameters);
population_proto_parameters_visualization(proto_parameters, Options);

close(figure(16), figure(17), figure(18));

%% Generating artificial cell lineages for one embryo and for the prototype

% with statistics from one embryo

% embryo number
i = 3;

% number of randomly generated cell lineages - warning: increasing the
% number of draws increases highly the computation time - ndraw = 300 gives
% good results
ndraw = 3;

mainmodelstat(i,ndraw, Data, Options, pop_parameters);

close(figure(19), figure(20), figure(21));

% with prototypical statistics
% number of randomly generated cell lineages - warning: increasing the
% number of draws increases highly the computation time - ndraw = 300 gives
% good results
ndraw = 3;
prototype(ndraw, Data, Options, proto_parameters);

close(figure(22), figure(23), figure(24));

clc