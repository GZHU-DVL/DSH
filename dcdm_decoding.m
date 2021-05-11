function dec_message = dcdm_decoding(lattice, Delta, dither, message, cosets, obs)

%Decoding of the message embedded in a sequence of watermarked signals using DC-DM
% Inputs:
%  lattice: type of shaping lattice: 'rep_coding', 'hexagonal', 'E8'
%  Delta: scaling factor for the shaping lattice
%  dither: secret dither
%  message: the message originally embedded
%  cosets: the minimum-norm vectors defining the cosets of the nested lattice code
%  obs: the matrix (n x No) with the watermarked vectors
%
% Outputs:
%  dec_message: the n-length decoded message


n = size(obs, 1);
No = size(obs, 2);
alphabet = [0 1];

dec_message = zeros(1,No);

switch lattice
    case 'rep_coding'
        
        for index_obs=1:No
            quant_error = zeros(n, size(cosets,2));
            for index_M=0:(size(cosets,2)-1)
                for index_dim=1:n
                    y = obs(index_dim,index_obs);
                    quant = Delta*round((y - Delta*index_M/size(cosets,2) - dither(index_dim))/Delta) + Delta*index_M/size(cosets,2) + dither(index_dim);
                    quant_error(index_dim, index_M+1) = quant - y;
                end
            end
            norm_error = [];
            for ind=1:size(cosets,2)
                norm_error(ind) = norm(quant_error(:,ind));
            end
            dec_message(index_obs) = find(norm_error==min(norm_error)) - 1;
        end
        
    otherwise   %'hexagonal' or 'E8' lattice
        
        for index_obs=1:No
            y = obs(:,index_obs);
            quant_error_norm = zeros(1,size(cosets,2));
            for index_coset = 1:size(cosets,2)
                quant = lattice_decoding((y' - cosets(:,index_coset)' - dither')/Delta , lattice)*Delta + ...
                cosets(:,index_coset)' + dither';
                quant_error_norm(index_coset) = norm(y - quant');
            end
            [min_error, index_min] = min(quant_error_norm);
            dec_message(index_obs) = index_min-1;
        end
        
end
