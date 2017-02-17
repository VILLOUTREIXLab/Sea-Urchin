function [ proto_parameters ] = prototypical_parameters( pop_parameters )
%PROTOTYPICAL_PARAMETERS computes prototypical statistics as centroid of
%the cohort using symmetrized Kullback-Leibler divergence

Ng = 12;% max number of generations
Nsp = 4;% number of cell populations

% intialization of the variables
meancyclelengthProto = zeros(Ng,Nsp,2);
meandivisiontimeProto = zeros(Ng,Nsp,2);
meanlogvolumeProto = zeros(Ng,Nsp,2);
meanlogcoefvolProto = zeros(Ng,Nsp,2);
meanlogsurfaceProto = zeros(Ng,Nsp,2);
meanlogcoefsurfProto = zeros(Ng,Nsp,2);


% loop over cell populations
for j =1:4,
    
    % loop over generations
    for g = 6:12,
        
        % division time
        mu = [];
        sigma = [];
        
        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meandivisiontime(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the division time mean values
            mu = [mu,pop_parameters.meandivisiontime(g,j,1,indtemp(i))];
            % aggregate the division time standard deviations
            sigma = [sigma,sqrt(pop_parameters.meandivisiontime(g,j,2,indtemp(i)))];
        end
        
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean division time
            meandivisiontimeProto(g,j,1) = val(1);
            % and prototypical division time standard deviation
            meandivisiontimeProto(g,j,2) = (val(2))^2;
        end
        
        % cell cycle length
        val = 0;
        mu = [];
        sigma = [];
        
        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meancyclelength(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the cell cycle length mean values
            mu = [mu,pop_parameters.meancyclelength(g,j,1,indtemp(i))];
            % aggregate the cell cycle length standard deviations
            sigma = [sigma,sqrt(pop_parameters.meancyclelength(g,j,2,indtemp(i)))];
        end
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean cell cycle length
            meancyclelengthProto(g,j,1) = val(1);
            % and prototypical cell cycle length standard deviation
            meancyclelengthProto(g,j,2) = (val(2))^2;
        end
        
        % log volume
        val = 0;
        mu = [];
        sigma = [];
        
        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meanlogvolume(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the log volume mean values
            mu = [mu,pop_parameters.meanlogvolume(g,j,1,indtemp(i))];
            % aggregate the log volume standard deviations
            sigma = [sigma,sqrt(pop_parameters.meanlogvolume(g,j,2,indtemp(i)))];
        end
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean log volume
            meanlogvolumeProto(g,j,1) = val(1);
            % and prototypical log volume standard deviation
            meanlogvolumeProto(g,j,2) = (val(2))^2;
        end
        
        % log daughter/mother volume ratio
        val = 0;
        mu = [];
        sigma = [];

        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meanlogcoefvol(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the log daughter/mother volume ratio mean values
            mu = [mu,pop_parameters.meanlogcoefvol(g,j,1,indtemp(i))];
            % and standard deviation
            sigma = [sigma,sqrt(pop_parameters.meanlogcoefvol(g,j,2,indtemp(i)))];
        end
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean log coef mother/daughter volume ratio
            meanlogcoefvolProto(g,j,1) = val(1);
            % and prototypical log coef mother/daughter volume ratio
            % standard deviation
            meanlogcoefvolProto(g,j,2) = (val(2))^2;
        end
        
        % log surface area
        val = 0;
        mu = [];
        sigma = [];
        
        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meanlogsurface(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the log surface area mean values
            mu = [mu,pop_parameters.meanlogsurface(g,j,1,indtemp(i))];
            % and standard deviation
            sigma = [sigma,sqrt(pop_parameters.meanlogsurface(g,j,2,indtemp(i)))];
        end
        
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean log surface area
            meanlogsurfaceProto(g,j,1) = val(1);
            % and prototypical log surface area standard deviation
            meanlogsurfaceProto(g,j,2) = (val(2))^2;
        end
        
        % log daughter/mother surface ratio
        val = 0;
        mu = [];
        sigma = [];
        
        % find the indices of the embryos having statistics for the current
        % population and generation
        indtemp = find(pop_parameters.meanlogcoefsurf(g,j,1,:));
        for i = 1:length(indtemp);
            % aggregate the log daughter/mother surface area ratio mean
            % values
            mu = [mu,pop_parameters.meanlogcoefsurf(g,j,1,indtemp(i))];
            % and standard deviation
            sigma = [sigma,sqrt(pop_parameters.meanlogcoefsurf(g,j,2,indtemp(i)))];
        end
        
        if length(indtemp)>0,
            % compute the centroid of the cohort using symmetrized Kullback
            % Leibler
            val = centroidBregman(mu,sigma);
            % prototypical mean log daughter/mother surface area ratio
            meanlogcoefsurfProto(g,j,1) = val(1);
            % and prototypical standard deviation
            meanlogcoefsurfProto(g,j,2) = (val(2))^2;
        end
    end
end

% storing the results
proto_parameters.meancyclelengthProto = meancyclelengthProto;
proto_parameters.meandivisiontimeProto = meandivisiontimeProto;
proto_parameters.meanlogvolumeProto = meanlogvolumeProto;
proto_parameters.meanlogcoefvolProto = meanlogcoefvolProto;
proto_parameters.meanlogsurfaceProto = meanlogsurfaceProto;
proto_parameters.meanlogcoefsurfProto = meanlogcoefsurfProto;

end

