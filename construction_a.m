function cosets = construction_a(p, g, G, debug)

%Constructs a nested lattice code using 'construction A'. The inputs are a scalar
%p (prime number?, if yes, ^ is a good lattice for MSE quantization) that indicates 
% the number of codewords, a generating vector g (size n x 1), a matrix 
%generator G of size n x n (where the basis vectors are given in rows) for the coarse
%lattice and a parameter (debug) which when set to other value than 0 plots
%the lattice points and the corresponding Voronoi regions (works only in 2
%dimensions).
%The output is a n x p matrix containing the cosets of the nested lattice code.

% clear all
% close all
% clc;
% %Example using a hexagonal lattice as coarse lattice:
% p = 9;
% % g = [1 2]';
% g = [2 3]'; %generating vector for the fine lattice
% Delta = 1;
% Ghex = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix
% G=Ghex;
% % cosets = construction_a(p, g, Ghex, 1)
% debug = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %generating vector
% g = rand(n,1)*p;
% g = floor(g);   %the elements of g are uniformly distributed in {0, 1, ..., p-1}
x=zeros(2, p);
for q=0:p-1
    x(:,q+1) = mod(q*g, p)/p;
end
cosets = x; %codewords or coset representatives (column vectors)
cosets*p;
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
%     figure(1), plot(x(1,:), x(2,:), '*'), grid on, title('Original lattice points generated by Construction A')
%     figure(2), voronoi(x(1,:)*p, x(2,:)*p), axis equal, axis([0, 1, 0, 1]*p), grid on,
%     title('Original lattice points (in the unit cube) generated by Construction A')
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
% figure(3), plot(y(:,1), y(:,2), '.'), hold on
% title('Fine lattice \Lambda_f (dots) and coarse lattice \Lambda_c (circles)')
% figure(3), hold on, plot(hex_latt(:,1), hex_latt(:,2), 'ro'), 
% hold off, axis equal
if (debug & length(g)==2)
   voronoi(y(:,1), y(:,2), 'r') % fine lattice points and lines
    [vx, vy] = voronoi(hex_latt(:,1), hex_latt(:,2));
    hold on, plot(hex_latt(:,1), hex_latt(:,2), '.k', 'MarkerSize', 18), % coarse lattice points
    plot(vx, vy, 'k', 'Linewidth', 2), hold off  % coarse lattice lines
    axis (1.5*[-1 1 -1 1]),
%     title('Fine lattice \Lambda_f (red) and coarse lattice \Lambda_c (black)')
%     axis equal
%     disp('Press any key ...')
%     pause
end

%this is for one of the cosets of the checkerboard lattice
%figure(8), plot(x=[-1.5 -0.5 0.5 1.5 -2.5 -1.5 -0.5 0.5 -0.5 0.5 1.5 2.5 0.5 1.5 2.5 3.5 -3.5 -2.5 -1.5 -0.5 -1.5 1.5], ...
% y=[-1.5 -0.5 0.5 1.5 -1.5 -0.5 0.5 1.5 -1.5 -0.5 0.5 1.5 -1.5 -0.5 0.5 1.5 -1.5 -0.5 0.5 1.5 1.5 -1.5], 'ok')
% voronoi(x,y), axis (2*[-1 1 -1 1])

