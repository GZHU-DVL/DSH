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
p = 9;
% g = [1 2]';
g = [2 3]'; %generating vector for the fine lattice
Delta = 1;
Ghex = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix
G=Ghex;
cosets = construction_a(p, g, Ghex, 0);
alpha=0.85;
% dither=rand_obs(1, Delta/3, lattice, dimension)'
dither=[0.1542 0.2657]';
% theta=2*pi*rand(1)
theta=0;
R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %generating vector
% g = rand(n,1)*p;
% g = floor(g);   %the elements of g are uniformly distributed in {0, 1, ..., p-1}
num_obs=30; % Monte Carlo number to approximate SER (symbol error rate)
No=num_obs;
debug = 2;
% hat_dither=zeros(2,Times);
obs =R*rand_obs(num_obs, Delta, lattice, dimension)';
message = floor(rand(1, No)*p); %sequence of symbols (p-ary, equiprobable) to be embedded
alphabet = 0:p-1;
message = mod(message - message(1), length(alphabet)); 
% watermark embedding
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
%dither exhaustive search begin
%namely, search precision(point) is  about 0.1 for Delta=1 lattice
% k=1;
point=2000;
% for point=100:100:1000 
vol_voronoi=(sqrt(3)*Delta.^2)/2;
vol_point=vol_voronoi/point;
step_search=sqrt(vol_point);
% four linear equations listed as follows from the 1- to 4-th quadrant
% respectively, y=-sqrt(3)*x+Delta; y=sqrt(3)*x+Delta; y=-sqrt(3)*x-Delta;
% y=sqrt(3)*x-Delta
pn=1;
est_dither=zeros(2,point);
%search divided into three regions
for x=-Delta/sqrt(3):step_search:(-Delta/(2*sqrt(3))-step_search)
    for y=(-sqrt(3)*x-Delta):step_search:(sqrt(3)*x+Delta)
        est_dither(:,pn)=[x, y]';
        % decoding using the estimated dither
        dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
        % compute the symbol error rate (SER)
       ser(pn)=sum(dec_message~=message)/No;
        pn=pn+1;
    end
end
for x=(-Delta/(2*sqrt(3))):step_search:(Delta/(2*sqrt(3)))
    for y=(-Delta/2):step_search:(Delta/2)
        est_dither(:,pn)=[x, y]';
        % decoding using the estimated dither
        dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
        % compute the symbol error rate (SER)
       ser(pn)=sum(dec_message~=message)/No;
        pn=pn+1;
    end
end
for x=(Delta/(2*sqrt(3))+step_search):step_search:(Delta/sqrt(3))
    for y=(sqrt(3)*x-Delta):step_search:(-sqrt(3)*x+Delta)
        est_dither(:,pn)=[x, y]';
        % decoding using the estimated dither
        dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
        % compute the symbol error rate (SER)
       ser(pn)=sum(dec_message~=message)/No;
        pn=pn+1;
    end
end
index=find(ser==0);
pro_dither=est_dither(:,index);
x=zeros(2, p);
for q=0:p-1
    x(:,q+1) = mod(q*g, p)/p;
end
x=cosets;
cosets = x; %codewords or coset representatives (column vectors)
cosets*p;
% % the below needs only to show the shape of lattice using construction_A
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
% axis equal
plot(dither(1,1), dither(2,1),'g.', 'MarkerSize', 18);
plot(pro_dither(1,:), pro_dither(2,:),'b.', 'MarkerSize', 4);
set(gca,'FontSize',18,'FontName','Times New Roman');
hold off;

