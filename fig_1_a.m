%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure showing a  nested-lattice code obtained through the self-similar construction 
% Hexagonal shaping lattice
% The rate of the code is log(9)/2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc;
lattice='hexagonal';
% num_obs=500;
Delta = 1.0954;  %scaling factor
G = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix for the hexagonal lattice
%% vol(v(^))=det(G)=Delta.^2*sqrt(3)/2
sqrt_order = sqrt(4); %the order of the partition is sqrt_order^2
debug = 1;
% obs =rand_obs(num_obs, Delta/3.3, lattice, 2);
[x, y] = meshgrid(-10:1:10);
combs = [];
for i=1:size(x,2)
    combs = [combs; [x(:,i), y(:,1)]];
end

hex_latt1 = combs*G/Delta;  %coarse lattice
hex_latt1=hex_latt1; % rotation by theta in the direction of the counterclockwise

[x, y] = meshgrid(-1:1:1);
combs = [];
for i=1:size(x,2)
    combs = [combs; [x(:,i), y(:,1)]];
end
coset_leaders = combs*G/Delta;  
coset_leaders = coset_leaders./sqrt_order;   %coset leaders of the nested lattice code

[vx, vy] = voronoi(hex_latt1(:,1), hex_latt1(:,2));
figure(1), plot(hex_latt1(:,1), hex_latt1(:,2), '.k', 'MarkerSize', 18),
hold on, plot(vx, vy, 'k', 'Linewidth', 1.5), hold off

hex_latt2 = hex_latt1/sqrt_order;   %fine lattice
[vx, vy] = voronoi(hex_latt2(:,1), hex_latt2(:,2));
%figure(1), hold on, plot(hex_latt2(:,1), hex_latt2(:,2), '.k', 'MarkerSize', 5),
hold on, plot(vx, vy, 'r', 'Linewidth', 0.5), hold off

coset = coset_leaders(1,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '+k', 'MarkerSize', 6), hold off

coset = coset_leaders(2,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '*k', 'MarkerSize', 6), hold off

coset = coset_leaders(3,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'ok', 'MarkerSize', 6), hold off

coset = coset_leaders(4,:);
hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
    hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'dk', 'MarkerSize', 6), hold off

% coset = coset_leaders(5,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '.k', 'MarkerSize', 1), hold off
% 
% 
% coset = coset_leaders(6,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '+k', 'MarkerSize', 7), hold off
% 
% coset = coset_leaders(7,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), '*k', 'MarkerSize', 6), hold off
% 
% coset = coset_leaders(8,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'xk', 'MarkerSize', 5), hold off
% 
% coset = coset_leaders(9,:);
% hold on, plot(hex_latt1(:,1) + repmat(coset(:,1), size(hex_latt1,1), 1), ...
%     hex_latt1(:,2) + repmat(coset(:,2), size(hex_latt1,1), 1), 'ok', 'MarkerSize', 6), hold off

axis equal
axis (1.5*[-1 1 -1 1]),
% axis ([-1.5 1.5 -1.5 1.5]);
set(gca,'FontSize',18,'FontName','Times New Roman');



