function [ val ] = centroidBregman( mu, sigma )
%CENTRDOIBREGMAN computes the Bregman centroid of a set of gaussians
%distributions using the symmetrized Kullback Leibler divergence

% the centroid is the gaussian distribution defined by its mean and
% standard deviation minimizing the symmetrized Kullback Leibler divergence
% to all distributions
if length(mu) == length(sigma),
    val = fminsearch(@(x) kullbackleiblerdistance(mu,sigma,x),[mean(mu),sqrt(mean(sigma.^2))]);
end

end

