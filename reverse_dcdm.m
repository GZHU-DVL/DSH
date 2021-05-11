function [est_host] = reverse_dcdm(lattice, Delta, alpha, dither, message, cosets, obs)

n = size(obs, 1);
No = size(obs, 2);

quant_error = zeros(n, No);
        
switch lattice
    case 'rep_coding'               
        
        for index_obs=1:No
            for index_dim=1:n
                y = obs(index_dim, index_obs);
                quant = Delta*round((y - message(index_obs)*Delta/2 - dither(index_dim))/Delta) + ...
                    message(index_obs)*Delta/2 + dither(index_dim);
                quant_error(index_dim, index_obs) = quant - y;
            end
        end
        est_host = obs - alpha*quant_error/(1-alpha);   %estimated host

    otherwise

        for index_obs=1:No
            y = obs(:,index_obs);
            %second, the quantization error with the decoded message is computed 
            quant = lattice_decoding((y' - cosets(:,message(index_obs)+1)' - dither')/Delta , lattice)*Delta + ...
                    cosets(:,message(index_obs)+1)' + dither';
            quant_error(:,index_obs) = quant' - y;
        end
        %now, the host is reversed using the decoded message
        est_host = obs - alpha*quant_error/(1-alpha);
        
end
