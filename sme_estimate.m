% function [volume_iter, dist_centroid, vector_obs, Delta] = sme_estimate(lattice, dimensions, alpha, num_obs, num_reps, original_centroid, debug, seed)
function [volume_iter, dist_centroid, vector_obs, Delta] = sme_estimate(lattice, dimensions, alpha, num_obs, num_reps, debug, seed)
%[volume_iter] = sme_estimate(lattice, dimensions, alpha, num_obs, num_reps, centroid, debug)
%
% Dither estimation in lattice-quatization data hiding using a
% set-membership algorithm (OVE, Optimal Volume Ellipsoid)
%
%Inputs: 
%   Lattice: the embedding lattice (cubic, hexagonal, Dn, Dndual, E7dual, E8)
%   dimensions: dimensionality of the embedding lattice
%   alpha: distortion compensation parameter
%   num_obs: number of observations
%   num_reps: number of recirculations of the data
%   centroid: location of the secret dither
%   debug: debugging plots (works only for the hexagonal lattice). Value between 0 and 2
%   seed: seed for the random number generation
%
%Outputs:
%   volume_iter: volume of the bounding ellipsoid
%   dist_centroid: Euclidean distance between the center of the i-th ellipsoid and the centroid
%   vector_obs: observation indexes (No)
%   Delta: scaling parameter of the lattice for ensuring G(\Lambda) = 1/12

% rand('state', seed)

%now, the scaling parameter and bounding hyperplanes (matrix Xk and vector gammak) 
%are calculated according to the embedding lattice
switch lattice
    case 'cubic'
        Delta = 1;
        vol = Delta^dimensions;
        lattice_type = 0;
        Xk = Delta*(1-alpha)*eye(dimensions);
        for index_gamma=1:size(Xk,2)
                gammak(index_gamma) = (norm(Xk(:,index_gamma))^2)/2;
        end
        cov_radius = Delta*sqrt(dimensions)/2  %cubic lattice
    case 'hexagonal'
        %OJO, usamos la matriz generadora Mhex = [0, Delta; Delta*sqrt(3)/2, Delta/2];
        dimensions = 2;
        Ghexagonal = 0.080188;
        M = [1 0; 1/2 sqrt(3)*1/2];
        vol = sqrt(det(M*M'));
        Delta = sqrt((1/12)/(Ghexagonal*vol));
        lattice_type = 1;
        Xk = Delta*(1-alpha)*[0, sqrt(3)/2, sqrt(3)/2; 1, 1/2, -1/2];
        for index_gamma=1:size(Xk,2)
                gammak(index_gamma) = (norm(Xk(:,index_gamma))^2)/2;
        end
        cov_radius = Delta/sqrt(3); %hexagonal lattice
    case 'sphere'
        Gsphere = (gamma(0.5*dimensions+1)^(2/dimensions)/((dimensions+2)*pi));
        volsphere = (pi^(dimensions/2)/gamma(dimensions/2 + 1))*(1/2)^(dimensions);  %volume of the sphere of radius 0.5;
        Delta = sqrt(volsphere^(-2/dimensions)/(Gsphere*12));
        cov_radius = Delta/2
        lattice_type = 2;
    case 'Dn'
        GDn = (1/12)*(dimensions==2) + 0.0078543*(dimensions==3) + 0.076603*(dimensions==4) + 0.075786*(dimensions==5);
        vol = 2;
        Delta = sqrt(vol^(-2/dimensions)/(GDn*12));
        cov_radius = Delta*(dimensions==2) + Delta*1*(dimensions==3) + ...
            (Delta*sqrt(dimensions/2)/sqrt(2))*(dimensions>3);
        lattice_type = 3;
    case 'Dndual'
        GDndual = 0.078543*(dimensions==3) + 0.076603*(dimensions==4) + 0.075625*(dimensions==5);
        vol = 0.5;
        Delta = sqrt(vol^(-2/dimensions)/(GDndual*12));
        Xk_cube = Delta*(1-alpha)*eye(dimensions);
        binary_matrix = dec2bin(0:(2^dimensions)/2 - 1, dimensions);
        for i=1:size(binary_matrix,1)
            for j=1:size(binary_matrix,2)
                Xk_octahedron(i,j) = str2num(binary_matrix(i,j));
            end
        end
        Xk_octahedron(find(Xk_octahedron==0)) = -1;
        Xk_octahedron = -(1-alpha)*Delta*0.5*Xk_octahedron';
%         if dimensions==3
%             Xk_octahedron = Delta*(1-alpha)*0.5*[1 1 1;
%                                                  1 1 -1;
%                                                  1 -1 1;
%                                                  -1 1 1]';
%         elseif dimensions==4
%             Xk_octahedron = Delta*(1-alpha)*0.5*[1 1 1 1;
%                                                  1 1 1 -1;
%                                                  1 1 -1 1;
%                                                  1 1 -1 -1
%                                                  1 -1 1 1;
%                                                  1 -1 1 -1;
%                                                  1 -1 -1 1;
%                                                  1 -1 -1 -1]';
%         end
        Xk = [Xk_cube, Xk_octahedron];
        for index_gamma=1:size(Xk,2)
                gammak(index_gamma) = (norm(Xk(:,index_gamma))^2)/2;
        end
        cov_radius = (Delta*sqrt(5/3)*sqrt(3)/4)*(dimensions==3) + ...
            (Delta*sqrt(dimensions)/(2*sqrt(2)))*((dimensions>=4)&&(mod(dimensions,2)==0)) + ...
            (Delta*sqrt(2*dimensions-1)/4)*((dimensions>=5)&&(mod(dimensions,2)~=0));
        lattice_type = 4;
    case 'E7dual'
        dimensions = 7;
        GE7dual = 0.073116;
        vol = sqrt(0.5);
        Delta = sqrt(vol^(-2/dimensions)/(GE7dual*12));
        cov_radius = sqrt(7/8)*Delta;
        lattice_type = 6;
    case 'E8'
        dimensions = 8;
        GE8 = 0.071682;
        vol = 1;
        Delta = sqrt(vol^(-2/dimensions)/(GE8*12));
        %Delta^dimensions*(1-alpha)^dimensions, pause   %volume of the scaled voronoi cell
        cov_radius = Delta*1;
        lattice_type = 5;
        %first, the matrix with the even number of coordinates
        binary_matrix = dec2bin(0:(2^dimensions)-1, dimensions);
        for i=1:size(binary_matrix,1)
            for j=1:size(binary_matrix,2)
                aux_matrix(i,j) = str2num(binary_matrix(i,j));
            end
        end
        aux_matrix(find(aux_matrix==0)) = -1;
        aux_matrix = -aux_matrix;
        Xk1 = [];
        for i=1:size(aux_matrix,1)
            num_minus = length(find(aux_matrix(i,:)==-1));
            if mod(num_minus,2)==0
                Xk1 = [Xk1; aux_matrix(i,:)];
            end
        end
        Xk1 = (1-alpha)*Delta*0.5*Xk1';
        Xk2 = [];
        %now, permutations of (\pm 1^2,0^6)
        aux_matrix = [ones(dimensions-1,1), eye(dimensions-1)];
        aux_matrix2 = [];
        for i=2:size(aux_matrix,2)  %the columns of aux_matrix are interchanged
            aux2 = zeros(dimensions-1,dimensions);
            aux2(:,1:i-1) = aux_matrix(:,dimensions-(i-2):dimensions);
            aux2(:,i:dimensions) = aux_matrix(:,1:dimensions-(i-1));
            aux_matrix2 = [aux_matrix2; aux2];
        end
        Xk2 = [Xk2; (1-alpha)*Delta*unique([aux_matrix; aux_matrix2], 'rows')];
        %now, permutations of (\pm 1,-1,0^6)
        aux_matrix = [ones(dimensions-1,1), -eye(dimensions-1)];
        aux_matrix2 = [];
        for i=2:size(aux_matrix,2)  %the columns of aux_matrix are interchanged
            aux2 = zeros(dimensions-1,dimensions);
            aux2(:,1:i-1) = aux_matrix(:,dimensions-(i-2):dimensions);
            aux2(:,i:dimensions) = aux_matrix(:,1:dimensions-(i-1));
            aux_matrix2 = [aux_matrix2; aux2];
        end
        Xk2 = [Xk2; (1-alpha)*Delta*unique([aux_matrix; aux_matrix2], 'rows')];
        Xk = [Xk1, Xk2'];
        for index_gamma=1:size(Xk,2)
                gammak(index_gamma) = (norm(Xk(:,index_gamma))^2)/2;
        end
end
%generation of the sequence of observations and secret dither vector T
%called as original_centoid
% disp('original_centroid:');
original_centroid=rand_obs(1, Delta, lattice, dimensions);
% test applicability
% indices = [2 3 4 9 10 11 12 17 18 19 25 26];
% observations =real_obs(num_obs, indices)+repmat(original_centroid, num_obs,1);

observations = rand_obs(num_obs, Delta*(1-alpha), lattice, dimensions) + repmat(original_centroid, num_obs,1);
observations = repmat(observations, num_reps, 1);
num_obs = num_reps*num_obs(end);
% %permutation of the components
% broza_prov = observations;
% observations(:,1) = broza_prov(:,2);
% observations(:,2) = broza_prov(:,1);
% Y=X+T; X=Y-T
offset = observations(1,:); %translation to the center of the Voronoi region
dec_centroid = lattice_decodingc((original_centroid - offset)/Delta, lattice_type);
centroid = (original_centroid - offset) - dec_centroid*Delta;
if debug==2
    figure(10), plot(observations(:,1), observations(:,2), '.'),
    hold on, plot(centroid(1), centroid(2), 'rx', 'markersize', 16, 'Linewidth', 3), hold off
    title('observations')
%     axis equal
%     pause
end
volsphere = (pi^(dimensions/2)/gamma(dimensions/2 + 1));    %volume of the sphere with unit radius
%initial ellipsoid. It is a sphere of radius "(1-alpha)*cov_radius" centered at the origin
P_iter = (1-alpha)^2*cov_radius^2*eye(dimensions); %matrix for the initial ellipsoid in each iteration
sigma2k_iter = 1;
theta_iter = zeros(dimensions,1); %center of the initial ellipsoid
error = [];
est_theta = [];
dist_theta = [];
detP = [];
if debug 
    ellipsoids_pointer = figure; 
    figure(ellipsoids_pointer),
    patch([Delta/sqrt(3), Delta/(2*sqrt(3)), -Delta/(2*sqrt(3)), -Delta/sqrt(3), -Delta/(2*sqrt(3)), Delta/(2*sqrt(3))], ...
        [0, Delta/2, Delta/2, 0, -Delta/2, -Delta/2], [0.7, 0.7, 0.7])
end
if debug==2 fig_pointer = figure;   end
volume_iter = zeros(1,size(observations, 1));   %volume of the bounding ellipsoids
volume_iter(1) = sqrt(det(P_iter))*volsphere;
dist_centroid = zeros(1,size(observations, 1)); %Euclidean distance to the centroid
dist_centroid(1) = norm(centroid - theta_iter');
for obs_index=2:size(observations, 1)
    obs_index;
    if debug
        [V, Lambda] = eig(P_iter*sigma2k_iter);
        aux = sum(Lambda,1);
        majoraxis = sqrt(aux(find(aux==max(aux),1)));
        minoraxis = sqrt(aux(find(aux==min(aux),1)));
        eccentricity = sqrt(1 - minoraxis^2/majoraxis^2);
        aux_axis = V(:,find(aux==max(aux),1));
        aux_axis = aux_axis*sign(aux_axis(1)+eps);
        az = sign(aux_axis(2))*acos(aux_axis'*[1;0]/(norm(aux_axis)))*180/pi;
        theta_dec = lattice_decodingc((theta_iter + offset')/Delta, lattice_type);
        theta_desp = (theta_iter + offset') - theta_dec*Delta;
        circle_obs = ellipse1(theta_desp(1), theta_desp(2), [majoraxis, eccentricity], az, [], [], 'degrees', 1000);
        figure(ellipsoids_pointer), hold on, plot(circle_obs(:,1), circle_obs(:,2)),
        plot(original_centroid(1), original_centroid(2), 'rx', 'markersize', 16, 'Linewidth', 3), hold off, axis equal
        title('sequence of bounding ellipsoids')
    end
    center = observations(obs_index,:) - observations(1,:);    %translation of the observation
    dec_center = lattice_decodingc(center/Delta, lattice_type);
    center_iter = center - dec_center*Delta;                        
            %the parameters of the current ellipsoid
            P = P_iter;
            theta = theta_iter;
            sigma2k = sigma2k_iter;
            %shifted and scaled voronoi region
            if debug==2
                x0 = center(1);
                y0 = center(2);
                figure(fig_pointer)
                patch(x0 + (1-alpha)*[Delta/sqrt(3), Delta/(2*sqrt(3)), -Delta/(2*sqrt(3)), -Delta/sqrt(3), -Delta/(2*sqrt(3)), Delta/(2*sqrt(3))], ...
                    y0 + (1-alpha)*[0, Delta/2, Delta/2, 0, -Delta/2, -Delta/2], [0.7 0.7 0.7])
                axis equal
                hold on,
                plot(linspace(-1, 1, 20), (gammak(1) + Xk(2,1)*y0)*ones(1,20)/Xk(2,1), 'r')
                plot(linspace(-1, 1, 20), (-gammak(1) + Xk(2,1)*y0)*ones(1,20)/Xk(2,1), 'r')
                x = linspace(-0.5,0.5,20);
                y = (gammak(2) - Xk(1,2)*x + Xk(1,2)*x0 + Xk(2,2)*y0)/Xk(2,2);
                plot(x, y, 'k')
                y = (-gammak(2) - Xk(1,2)*x + Xk(1,2)*x0 + Xk(2,2)*y0)/Xk(2,2);
                plot(x, y, 'k')                
                y = (gammak(3) - Xk(1,3)*x + Xk(1,3)*x0 + Xk(2,3)*y0)/Xk(2,3);
                plot(x, y, 'g')
                y = (-gammak(3) - Xk(1,3)*x + Xk(1,3)*x0 + Xk(2,3)*y0)/Xk(2,3);
                plot(x, y, 'g')                
                plot(circle_obs(:,1) - offset(1), circle_obs(:,2) - offset(2)),
                hold off
                title(sprintf('observation %d\n', obs_index))
%                 pause
            end           
            %RECURSION            
            yk = center_iter*Xk;            
            color_ellipsoid = {'r', 'k', 'g'};
            for index_order = 1:size(Xk,2)
                %index_order
                [P, theta, sigma2k] = bound_OVE(Xk(:,index_order), yk(index_order), gammak(index_order), P, theta, sigma2k);
                %we check wether the centroid is inside the new ellipsoid or not
                if abs((centroid' - theta)'*inv(P)*(centroid' - theta)) > 1
                    fprintf(1, 'Error: centroid not in bounding ellipsoid')
                    %return
                end
                %new ellipsoid
                if debug==2
                    [V, Lambda] = eig(P*sigma2k);
                    aux = sum(Lambda,1);
                    majoraxis = sqrt(aux(find(aux==max(aux),1)));
                    minoraxis = sqrt(aux(find(aux==min(aux),1)));
                    eccentricity = sqrt(1 - minoraxis^2/majoraxis^2);
                    aux_axis = V(:,find(aux==max(aux),1));
                    aux_axis = aux_axis*sign(aux_axis(1)+eps);
                    az = sign(aux_axis(2))*acos(aux_axis'*[1;0]/(norm(aux_axis)))*180/pi;
                    circle = ellipse1(theta(1), theta(2), [majoraxis, eccentricity], az, [], [], 'degrees', 1000);
                    figure(fig_pointer), hold on, plot(circle(:,1), circle(:,2), color_ellipsoid{index_order}),
                    plot(centroid(1), centroid(2), 'rx', 'markersize', 12, 'Linewidth', 2), hold off
%                     pause
                end
            end
            if debug==2 close(fig_pointer); end                                
            %the computed ellipsoid is stored
            P_iter = P;
            theta_iter = theta;
            sigma2k_iter = sigma2k;
            volume_iter(obs_index) = sqrt(det(P_iter))*volsphere;
            dist_centroid(obs_index) = norm(centroid - theta_iter');
            %volume_iter(obs_index)                    
end %observations
vector_obs = 1:num_obs;
%volume_iter = [Delta^dimensions*vol, volume_iter];
% return
% %%%plots: volume of the bounding intersection and Euclidean distance to the centroid
% figure, semilogy(vector_obs, volume_iter)
% title(sprintf('%s lattice: %d dimensions, %d observations, %d recirculations\n', lattice, dimensions, num_obs/num_reps, num_reps))
% xlabel('observation index (N_o)'), ylabel('volume')
% figure, plot(vector_obs(2:end), dist_centroid)
% xlabel('observation index (N_o)'), ylabel('distance')
% title(sprintf('%s lattice: %d dimensions, %d observations, %d recirculations\n', lattice, dimensions, num_obs/num_reps, num_reps))



