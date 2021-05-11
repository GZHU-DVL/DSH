%This script launches simulations for testing SME and minimization algorithms for dither estimation
%It works for the Known Message Attack case
%UNCOMMENT THIS FOR USING IN LINUX
% %Sedumi is added to the path
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad/LMI/SeDuMi_1_1');
% %Yalmip is added to the path
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad/LMI/yalmip');
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad/LMI/yalmip/extras');
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad/LMI/yalmip/solvers');
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad/LMI/yalmip/demos');
% path(path, '/home/gts/alumnos/lpfreire/watermarking/seguridad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all
clc;
% lattice = 'E8'  %embedding lattice
lattice = 'hexagonal';
dimensions = 2; %dimensionality of the dither vector
num_obs = [5:5:50];%,15:5:30,40,50];  %number of observations
num_reps = 1;   %number of recirculations of the observed data
centroid = zeros(1,dimensions);  %the true centroid
alpha = 0.7   %distortion compensation parameter
debug = 0  %indicates whether debugging information must be plotted at the end or not (only works for the hexagonal lattice)
% method = 'inner_polytope'   %method for computing the bound to the intersection region
method = 'OVE'
%method = 'inner_ellipsoid'
num_realizations = 1000;   %number of realizations for averaging the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volume_iter = zeros(num_realizations, length(num_obs));
dist_centroid = zeros(num_realizations, length(num_obs));
for i=1:num_realizations
     switch method
        case 'OVE'
            for index=1:length(num_obs)
                j = num_obs(index);
                [volume, dist, vector_obs, Delta] = sme_estimate(lattice, dimensions, alpha, j, num_reps, debug, i-1);
                volume_iter(i,index) = volume(end);
                dist_centroid(i,index) = dist(end);
            end
        case 'inner_ellipsoid'
            [volume, vector_obs, Delta] = inner_ellipsoid(lattice, dimensions, alpha, num_obs, centroid, debug, i-1, 'inner');
            volume_iter(i,:) = volume;
        case 'inner_polytope'
            [volume, dist, vector_obs, Delta] = inner_polytope(lattice, dimensions, alpha, num_obs, centroid, debug, i-1);
            volume_iter(i,:) = volume;
            dist_centroid(i,:) = dist;
    end
  end
mean_dist = mean(dist_centroid.^2, 1);
%     Consumed=cputime-TS
figure
semilogy(num_obs/10, mean_dist, '-.k', 'Linewidth', 1);
legend('DSH','FontSize',18,'FontName','Times New Roman');
xlabel('Consumed Time (Seconds)','FontSize',18,'FontName','Times New Roman');
ylabel('MSE','FontSize',18,'FontName','Times New Roman')
set(gca,'FontSize',18,'FontName','Times New Roman');


