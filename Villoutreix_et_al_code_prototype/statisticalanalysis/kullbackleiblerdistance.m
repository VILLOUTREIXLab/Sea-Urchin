function [ res ] = kullbackleiblerdistance( mu,sigma, Proto )
%KULLBACKLEIBLERDISTANCE computes the symmetrized Kullback-Leibler
%divergence between a set of distributions and a distribution defined as
%Proto - it is assumed that all distributions are gaussian

% mean and standard deviation of the distribution Proto
muP = Proto(1);
sigmaP = Proto(2);

res = 0;
% for each distribution
for i = 1:length(mu),  
    % symmetrized Kullback Leibler divergence between the distribution and
    % Proto
    temp = 0.25*((sigma(i)^2 / sigmaP^2) + (sigmaP^2 / sigma(i)^2) + (muP - mu(i))^2*(1/sigma(i)^2 + 1/sigmaP^2) - 2);
    res = res + temp;
end

if (length(mu)>0),
    % the results are averaged
   res = res/length(mu); 
end

end

