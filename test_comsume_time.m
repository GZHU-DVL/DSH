%This script launches simulations for testing SME and minimization algorithms for dither estimation
%It works for the Known Message Attack case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
TS=cputime;
lattice = 'E8';  %embedding lattice
% lattice = 'hexagonal';
dimensions = 8; %dimensionality of the dither vector
% dimensions = 2; %dimensionality of the dither vector
% num_obs = [2:2:10,15,25,30,40,50];  %number of observations
num_obs = 700;
num_reps = 1;   %number of recirculations of the observed data
centroid = zeros(1,dimensions);  %the true centroid
alpha = 0.7;   %distortion compensation parameter
debug = 0;  %indicates whether debugging information must be plotted at the end or not (only works for the hexagonal lattice)
method = 'OVE';
num_realizations =5;   %number of realizations for averaging the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_iter = zeros(num_realizations, length(num_obs));
dist_centroid = zeros(num_realizations, length(num_obs));
for i=1:num_realizations
                [volume, dist, vector_obs, Delta] = sme_estimate(lattice, dimensions, alpha,num_obs, num_reps, debug, i-1);
                volume_iter(i) = volume(end);
                dist_centroid(i) = dist(end);    
end
mean_dist = mean(dist_centroid.^2, 1)
Consumed=cputime-TS


