function obs = rand_obs(num_obs, Delta, lattice, dimension)

%generation of "num_obs" random vectors uniformly distributed in the Voronoi
%region of the lattice "lattice" of dimensionality "dimension", scaled by
%"Delta"

switch lattice
    case 'cubic'
        obs = Delta*rand(num_obs,dimension) - Delta/2;
        
    case 'hexagonal'    %the hexagonal lattice in 2 dimensions
        cov_radius = Delta/sqrt(3);
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, 2);

        for counter=1:num_obs
            a = rand(1)*5000 - 2500;
            b = rand(1)*5000 - 2500;
           decoded_point = lattice_decodingc([a, b]/Delta, 1);
           %decoded_point = lattice_decoding([a, b]/Delta, 1);
            decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = [a b] - decoded_point;
        end
%     plot(obs(:,1),obs(:,2), '*');
    case 'sphere'
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);
        while counter<num_obs
            counter2 = counter2 + 1;
            aux = Delta*rand(1,dimension) - Delta/2;
            decoded_point = lattice_decodingc(aux/Delta, 2);
            %decoded_point = decoded_point*Delta;    %this is the quantized value
            if (sum(abs(decoded_point))==0)
                counter = counter + 1;
                obs(counter,:) = aux;
            end
        end
        
    case 'Dn'
        cov_radius = Delta*(dimension==2) + Delta*(1*(dimension==3) + Delta*(sqrt(dimension/2)/sqrt(2))*(dimension>3));
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);
        
        for counter=1:num_obs
            rand_point = rand(1,dimension)*5000 - 2500;
            %decoded_point = lattice_decoding(rand_point/Delta, 'Dn');
            decoded_point = lattice_decodingc(rand_point/Delta, 3);
            decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = rand_point - decoded_point;
        end


    case 'Dndual'
        cov_radius = (Delta*sqrt(5/3)*sqrt(3)/4)*(dimension==3) + ...
            (Delta*sqrt(dimension)/(2*sqrt(2)))*((dimension>=4)&&(mod(dimension,2)==0)) + ...
            (Delta*sqrt(2*dimension-1)/4)*((dimension>=5)&&(mod(dimension,2)~=0));
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);

        for counter=1:num_obs
            rand_point = rand(1,dimension)*5000 - 2500;
            %decoded_point = lattice_decoding(rand_point/Delta, 'Dndual');
            decoded_point = lattice_decodingc(rand_point/Delta, 4);
            decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = rand_point - decoded_point;
        end
         

    case 'An'
        a = floor(dimension+1)/2;
        cov_radius = sqrt(2*a*(dimension+1-a)/(dimension+1))/sqrt(2)
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);
        if dimension==1
            M = [1 -1]; 
        else
            M = [diag(ones(1,dimension),0) + diag(-1*ones(1,dimension-1),1), ...
                [zeros(dimension-1); -1]];
        end
        %M
        %M = [1 0 -1; 1/sqrt(3) -2/sqrt(3) 1/sqrt(3)]

        for counter=1:num_obs
            rand_point = rand(1,dimension)*5000 - 2500;
            trans_point = rand_point*M;
            decoded_point = lattice_decoding(trans_point/Delta, 'An');
            %decoded_point = lattice_decodingc(rand_point/Delta, 5);
            decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = rand_point - decoded_point*M'/2;
        end


    case 'E7dual'
        cov_radius = sqrt(7/8);
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);
        
        cosets = [0 0 0 0 0 0 0;
            0 0 0 0 1 1 1;
            0 0 1 1 0 0 1;
            0 0 1 1 1 1 0;
            0 1 0 1 0 1 0;
            0 1 0 1 1 0 1;
            0 1 1 0 0 1 1;
            0 1 1 0 1 0 0;
            1 0 0 1 0 1 1;
            1 0 0 1 1 0 0;
            1 0 1 0 0 1 0;
            1 0 1 0 1 0 1;
            1 1 0 0 0 0 1;
            1 1 0 0 1 1 0;
            1 1 1 1 0 0 0;
            1 1 1 1 1 1 1];
        
        for counter=1:num_obs
            rand_point = rand(1,dimension)*5000 - 2500;
            %decoded_point = lattice_decoding(rand_point/Delta, 'E7dual');
            decoded_point = lattice_decodingc(rand_point/Delta, 6, cosets');
            decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = rand_point - decoded_point;
        end
        
        
    case 'E8'
        cov_radius = 1;
        counter = 0;
        counter2 = 0;
        obs = zeros(num_obs, dimension);

        for counter=1:num_obs
            rand_point = rand(1,dimension)*20000 - 10000;
            decoded_point = lattice_decoding(rand_point/Delta, 'E8');
            %decoded_point = lattice_decoding(rand_point/Delta, 5);
           decoded_point = decoded_point*Delta;    %this is the quantized value
            obs(counter,:) = rand_point - decoded_point;
        end

end

