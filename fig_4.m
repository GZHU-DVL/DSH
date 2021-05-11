%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure showing a  nested-lattice code rotated by theta obtained through the self-similar construction 
% Hexagonal shaping lattice
% The rate of the code is log(9)/2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc;
lattice='hexagonal';
num_obs=500;
Delta = 1.0954;  %scaling factor
G = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix for the hexagonal lattice
%% vol(v(^))=det(G)=Delta.^2*sqrt(3)/2
sqrt_order = sqrt(9); %the order of the partition is sqrt_order^2
% theta=rand(1)*2*pi;
theta=0;
debug = 1;
% rotation matrix for two dimensions, det(R)=1 
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
obs =rand_obs(num_obs, Delta/3.3, lattice, 2);
[x, y] = meshgrid(-10:1:10);
combs = [];
for i=1:size(x,2)
    combs = [combs; [x(:,i), y(:,1)]];
end

hex_latt1 = combs*G/Delta;  %coarse lattice
hex_latt1=hex_latt1*R'; % rotation by theta in the direction of the counterclockwise

[x, y] = meshgrid(-1:1:1);
combs = [];
for i=1:size(x,2)
    combs = [combs; [x(:,i), y(:,1)]];
end

coset_leaders = combs*G/Delta;  
coset_leaders = coset_leaders*R'./sqrt_order;   %coset leaders of the nested lattice code

[vx, vy] = voronoi(hex_latt1(:,1), hex_latt1(:,2));
figure(1), plot(hex_latt1(:,1), hex_latt1(:,2), '.k', 'MarkerSize', 18),
hold on, plot(vx, vy, 'k', 'Linewidth', 2), hold off
axis equal
axis ([-0.6 .6 -.6 0.6]);
hex_latt2 = hex_latt1/sqrt_order;   %fine lattice
[vx, vy] = voronoi(hex_latt2(:,1), hex_latt2(:,2));
hold on, plot(vx, vy, 'r', 'Linewidth', 0.5), hold off

% coset = coset_leaders(1,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'sk', 'MarkerSize', 6), hold off
% 
% coset = coset_leaders(2,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'hk', 'MarkerSize', 6), hold off
% 
% coset = coset_leaders(3,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'vk', 'MarkerSize', 6), hold off

coset = coset_leaders(4,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'dk', 'MarkerSize', 6), hold off

% the following line added by wyg
hold on, plot(repmat(coset(:,1), 500, 1)+obs(:,1), ...
    repmat(coset(:,2), 500, 1)+obs(:,2), '*g', 'MarkerSize', 2), hold off

% coset = coset_leaders(5,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '.k', 'MarkerSize', 1), hold off
% 
% 
% coset = coset_leaders(6,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '+k', 'MarkerSize', 7), hold off

coset = coset_leaders(7,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '*k', 'MarkerSize', 6), hold off
% the following line added by wyg
hold on, plot(repmat(coset(:,1), 500, 1)+obs(:,1), ...
    repmat(coset(:,2), 500, 1)+obs(:,2), '*b', 'MarkerSize', 2), hold off

% coset = coset_leaders(8,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'xk', 'MarkerSize', 5), hold off
% 
% coset = coset_leaders(9,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'ok', 'MarkerSize', 6), hold off

axis equal
axis (0.6*[-1 1 -1 1]),
set(gca,'FontSize',18,'FontName','Times New Roman');

