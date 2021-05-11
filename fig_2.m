% test the impact of alpha in DC-DM by Constructs a nested lattice code using 'construction A'. The inputs are a scalar
%p (prime number?, if yes, ^ is a good lattice for MSE quantization) that indicates 
% the number of codewords, a generating vector g (size n x 1), a matrix 
%generator G of size n x n (where the basis vectors are given in rows) for the coarse
%lattice and a parameter (debug) which when set to other value than 0 plots
%the lattice points and the corresponding Voronoi regions (works only in 2
%dimensions).

clear all
close all
clc;
%Example using a hexagonal lattice as coarse lattice:
lattice= 'hexagonal';
lattice_type=1;
dimension=2;
num_obs=5000;
No=num_obs;
p = 9;
% g = [1 2]';
g = [2 3]'; %generating vector for the fine lattice
Delta = 1;
Ghex = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix
G=Ghex;
cosets = construction_a(p, g, Ghex, 0);
alpha=1;
dither=rand_obs(1, Delta, lattice, dimension)'
% dither=[0 0]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 2;
obs =rand_obs(num_obs, Delta, lattice, dimension)';
message = floor(rand(1, No)*p); %sequence of symbols (p-ary, equiprobable) to be embedded
alphabet = 0:p-1;
message = mod(message - message(1), length(alphabet)); 
for index_obs = 1:No
    host = obs(:,index_obs);
    coset_index = message(index_obs)+1;    
%     quant = lattice_decodingc((host' - cosets(:,coset_index)' - dither(:,index_obs)')/Delta , lattice_type)*Delta + ...
%         cosets(:,coset_index)' + dither(:,index_obs)'; % dither(:,index_obs)' take place of dither by wyg
    quant = lattice_decodingc((host' - cosets(:,coset_index)' - dither')/Delta , lattice_type)*Delta + ...
        cosets(:,coset_index)' + dither'; % dither(:,index_obs)' take place of dither by wyg
    watermark(:, index_obs) = alpha*(quant' - host);
    obs(:,index_obs) = host + watermark(:, index_obs);
end

x=zeros(2, p);
for q=0:p-1
    x(:,q+1) = mod(q*g, p)/p;
end

cosets = x; %codewords or coset representatives (column vectors)

% the below needs only to show the shape of lattice using construction_A
if (debug & length(g)==2)
    x = [x, x + repmat([1; 0], 1, p), x + repmat([0; 1], 1, p), x + repmat([-1; 0], 1, p), x + repmat([0; -1], 1, p), ...
        x + repmat([1; -1], 1, p), x + repmat([-1; 1], 1, p), x + repmat([1; 1], 1, p), x + repmat([-1; -1], 1, p), ...
        x + repmat([0; -2], 1, p), x + repmat([0; 2], 1, p), x + repmat([2; -2], 1, p), x + repmat([2; -1], 1, p), ...
        x + repmat([2; 0], 1, p), x + repmat([2; 1], 1, p), x + repmat([2; 2], 1, p), x + repmat([1; -2], 1, p), ...
        x + repmat([1; 2], 1, p), x + repmat([-2; -2], 1, p), x + repmat([-2; -1], 1, p), x + repmat([-2; 0], 1, p), ...
        x + repmat([-2; 1], 1, p), x + repmat([-2; 2], 1, p), x + repmat([-1; -2], 1, p), x + repmat([-1; 2], 1, p), ...
        x + repmat([-3; -3], 1, p), x + repmat([-3; -2], 1, p), x + repmat([-3; -1], 1, p), x + repmat([-3; 0], 1, p), ...
        x + repmat([-3; 1], 1, p), x + repmat([-3; 2], 1, p), x + repmat([-3; 3], 1, p), x + repmat([2; -3], 1, p)];
end
cosets = (cosets'*G)';   
% the below needs only to show the shape of lattice using construction_A
if (debug & length(g)==2)
    y = x'*G;    %points of the fine lattice
    combs = [0 0; 0 1; 0 -1; 1 0; 1 1; 1 -1; -1 0; -1 1; -1 -1; -2 -2; -2 -1; -2 0; -2 1; -2 2; 2 -2; 2 -1; 2 0; 2 1; 2 2; ...
        -3 -3; -3 -2; -3 -1; -3 0; -3 1; -3 2; -3 3; 3 -3; 3 -2; 3 -1; 3 0; 3 1; 3 2; 3 3; -1 -3; -1 -2; -1 2; -1 3; ...
        1 -3; 1 -2; 1 2; 1 3; -2 -3; 2 3; 0 -2; 0 2; 0 -3; 0 3];
    hex_latt = combs*G;
    hex_latt = hex_latt + repmat(cosets(:,1)', size(hex_latt, 1), 1);   %lattice_points as row vectors
end
if (debug & length(g)==2)
    figure(5), voronoi(y(:,1), y(:,2), 'r') % fine lattice points and lines
    [vx, vy] = voronoi(hex_latt(:,1), hex_latt(:,2));
    hold on, plot(hex_latt(:,1), hex_latt(:,2), '.k', 'MarkerSize', 18), % coarse lattice points
    plot(vx, vy, 'k', 'Linewidth', 2), hold off  % coarse lattice lines
    axis (0.6*[-1 1 -1 1]), 
end
hold on;
plot(obs(1,:), obs(2,:), 'b.');
axis equal
axis (0.6*[-1 1 -1 1]),
set(gca,'FontSize',18,'FontName','Times New Roman');
hold off;
