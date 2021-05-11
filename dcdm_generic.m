function [obs, host, cosets, alphabet, message, watermark, Delta, alpha, cov_radius, dither, DWR] = dcdm_generic(No_x, Delta, alpha, lattice, n, p, g, state, DWR)

%DC-DM with repetition coding of rate n. All the host vectors are
%watermarked using the same secret dither signal.
% 
%Inputs: number of observations (No_x), quantization step (Delta), distortion
%  compensation parameter (alpha), the type of shaping lattice (lattice), 
%  the dimensionality of the lattice (n), the order of the lattice partition 
%  or cardinality of the alphabet (p), the generating vector (g), 
%  the seed of the pseudorandom generator (state), and the DWR in dB (DWR).
%
%Outputs: a matrix n x No_x with the watermarked vectors (obs), 
%  the cosets of the nested lattice code (cosets), the elements
%  of the alphabet (alphabet) the sequence of
%  embedded symbols (message), the watermark vectors (watermark), the 
%  quantization step (Delta), the parameter alpha, the covering 
%  radius of the coarse lattice (cov_radius),
%  the secret dither (dither), and the DWR in dB (DWR).
%
%The lattices supported are 'hexagonal', 'E_8'
%input: (No_x, Delta, alpha, lattice, n, p, g, state, DWR)
%Usage example: [obs, cosets, alphabet, message, watermark, Delta, alpha, dither, DWR] = ...
%       dcdm_generic(1000, 1, 0.8, 'hexagonal', 1434,9,[1 2]',0,30)
%  the above is revised by WYG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%initialization of the pseudorandom generators
rand('state', state);
randn('state', state);

sigmax = 10;

%parameters for the hexagonal lattice
%n = 2;
%p = 9;
%g = [1 2]'; %generating vector for the fine lattice

switch lattice
    case 'hexagonal'
        n = 2;
        M = [0, 1; sqrt(3)/2, 1/2];% an arbitrary n*n real-valued matrix M defines a lattice ^, M also known as generating matrix
        vol = sqrt(det(M*M'));
        Ghexagonal = 0.080188;  %normalized second order moment
        %Delta = sqrt((1/12)/(Ghexagonal*vol))
        Delta = sqrt((sigmax^2*10^(-DWR/10)/(alpha^2))/(Ghexagonal*vol));
        G = [0, Delta; Delta*sqrt(3)/2, Delta/2];   %generator matrix of the hexagonal lattice (basis vectors in rows)        
        cov_radius = Delta/sqrt(3);
        lattice_type = 1;
        
%        Ghexagonal*vol^(2/n)*Delta^2    %second order moment per dimension of the coarse lattice
%        (1-alpha)^2*Ghexagonal*vol^(2/n)*Delta^2
        
    case 'E8'
        n = 8;
        GE8 = 0.071682; %normalized second order moment
        G = diag(ones(1,8), 0) + diag(-ones(1,7), -1);
        G(1,1) = 2;
        G(8,:) = 1/2;   %generator matrix of the lattice E_8 (basis vectors in rows)
        vol = 1;
        Delta = sqrt((sigmax^2*10^(-DWR/10)/(alpha^2))/(GE8*vol))
        G = Delta*G;
        cov_radius = Delta;
        %vol = Delta^n;
        lattice_type = 5;

%        GE8*vol^(2/n)*Delta^2    %second order moment per dimension of the coarse lattice
%        (1-alpha)^2*GE8*vol^(2/n)*Delta^2

        
end

cosets = construction_a(p, g, G, 0); %coset leaders of the fine lattice (we still need to compute their minimum norm versions)

for i=1:size(cosets,2)  
    quant = lattice_decodingc(cosets(:,i)'/Delta, lattice_type)*Delta;
    cosets(:,i) = cosets(:,i) - quant'; %minimum norm coset representatives (column vectors)
end
% plot(cosets(1,:), cosets(2,:), 'ro');  % added by WYG to plot cosets' distribution
if length(No_x)==1    %the first argument is the number of observations to be generated
    dither = rand_obs(No_x, Delta, lattice, n)'; 
else    %the first argument is the set of host vectors
    %dither = rand_obs(size(No_x,2), Delta, lattice, n)'; 
    dither = rand_obs(1, Delta, lattice, n)'; 
end
%generation of the secret dither, size(No_x,2)  is taken place of 1 by wyg
%dither = zeros(n,1);    %this simplifies computation of the distance between the dither and its estimate and does not imply any loss of generality
% the above line is removed by wyg
% plot(dither(1,:), dither(2,:), '*') % added by WYG to plot the random dither distribution
if length(No_x)==1    %the first argument is the number of observations to be generated
    No = No_x;
    x = randn(n, No)*sigmax;    %random generation of the host vectors
else    %the first argument is the set of host vectors
    No = size(No_x,2);
    x = No_x;   %the host vectors
end

message = floor(rand(1, No)*p); %sequence of symbols (p-ary, equiprobable) to be embedded
alphabet = 0:p-1;
%message = zeros(1, No);
message = mod(message - message(1), length(alphabet));  %the message corresponding to the first observation is set to 0
obs = zeros(n, No); %the watermarked vectors
watermark = zeros(n, No);   %the watermark
for index_obs = 1:No
    host = x(:,index_obs);
    coset_index = message(index_obs)+1;    
%     quant = lattice_decodingc((host' - cosets(:,coset_index)' - dither(:,index_obs)')/Delta , lattice_type)*Delta + ...
%         cosets(:,coset_index)' + dither(:,index_obs)'; % dither(:,index_obs)' take place of dither by wyg
    quant = lattice_decodingc((host' - cosets(:,coset_index)' - dither')/Delta , lattice_type)*Delta + ...
        cosets(:,coset_index)' + dither'; % dither(:,index_obs)' take place of dither by wyg
    watermark(:, index_obs) = alpha*(quant' - host);
    obs(:,index_obs) = host + watermark(:, index_obs);
end

mean_watermark = mean(watermark, 2);
cov_matrix_watermark = (watermark - repmat(mean_watermark, 1, No))*(watermark - repmat(mean_watermark, 1, No))'/No;
sigmaw2 = trace(cov_matrix_watermark)/n;    %mean embedding distortion per dimension
sigmax2 = mean(var(x, 1, 2));
DWR = 10*log10(sigmax2/sigmaw2);


host = x;
