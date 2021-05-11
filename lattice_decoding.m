function dec_point = lattice_decoding(x, lattice)
%this function decodes (quantizes) the point x according to the lattice
%"lattice". Lattices supported are Z^n, hexagonal, D_n, D_n dual, A_n, E_7 dual, E_8 

switch lattice
    case 'cubic'
        for i=1:length(x)
            [dec_point(i), w(i)] = rounding(x(i));
        end
        
        
    case 'hexagonal'
       %the hexagonal lattice can be constructed through the union of two rectangular lattices
        aux(1,1) = rounding(x(1)*sqrt(3)/3)*3/sqrt(3);
        aux(1,2) = rounding(x(2));
        offset = [3/sqrt(3), 1]/2;
        aux(2,1) = rounding((x(1)- offset(1))*sqrt(3)/3)*3/sqrt(3) + offset(1);
        aux(2,2) = rounding((x(2)- offset(2))) + offset(2);
        dists = [norm(x - aux(1,:)),  norm(x - aux(2,:))];
        index = find(dists==min(dists));
        if length(index)>1  %in case of a tie, the vector with minimum norm is chosen
            norm_decoded = sum(aux.^2, 2);
            dec_point = aux(find(norm_decoded==min(norm_decoded)),:);
        else
            dec_point = aux(index, :);
        end
        
        
    case 'sphere'
        for i=1:length(x)
            [dec_point(i), w(i)] = rounding(x(i));
        end
        if norm(x-dec_point)^2 > 1/4    %in this case, the point does not belong to any sphere
            dec_point = NaN;
        end
        
        
    case 'Dn'   %the checkerboard lattice (n>=3)
        for i=1:length(x)
            [f(i), w(i)] = rounding(x(i));
        end        
        aux = abs(f - x);
        index = find(aux==max(aux));
        g = f;
        g(index(end)) = w(index(end));
        if mod(sum(f),2)==0 %the lattice point with even sum is chosen
            dec_point = f;
        else
            dec_point = g;
        end
        
        
    case 'Dndual'  %the dual of the lattice D3. Dn* can be seen as the union of 2 cosets of Zn
        cosets = [zeros(1, length(x)); ones(1, length(x))/2];
        dist = zeros(1,2);
        for i=1:size(cosets,1)
            %aux(i,:) = lattice_decoding(x - cosets(i,:), 'cubic');
            aux(i,:) = lattice_decodingc(x - cosets(i,:), 0);
            aux(i,:) = aux(i,:) + cosets(i,:);
            dist(i) = norm(x - aux(i,:))^2;
        end
        index = find(dist==min(dist));
        dec_point = aux(index,:);
                                 
        
    case 'An'   %x \in \mathbb{R}^8
        M = [1 0 -1; 1/sqrt(3) -2/sqrt(3) 1/sqrt(3)];
        %M = [1 -1 0; 0 1 -1];
%        x = x*M;
        
        s = sum(x);
        x2 = x - s/length(x);
        for i=1:length(x)
            fx2(i) = rounding(x2(i));
        end
        deff = sum(fx2)    %deficiency
        dec_point = fx2;
        if deff==0
            %dec_point = dec_point*M'/2;
            return
        else
            delta = abs(deff-fx2);  %quantization error per dimension
            [sorted, indexes] = sort(delta);
            sortedfx2 = fx2(indexes);
            if deff>0
                dec_point(indexes(1:deff)) = fx2(indexes(1:deff)) - 1;
            else
                dec_point(indexes(end:-1:end+deff+1)) = fx2(indexes(end:-1:end+deff+1)) + 1;
            end
            %dec_point = dec_point*M'/2;
        end        
        
        
    case 'E7dual'   %E_7^* is the union of 16 cosets of 2Z^7. The cosets are given by the codewords of H(7,4)
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

        dist = zeros(1,size(cosets,1));
        for i=1:size(cosets,1)
            %aux(i,:) = lattice_decoding((x/sqrt(0.5) - cosets(i,:))/2, 'cubic');
            aux(i,:) = lattice_decodingc((x/sqrt(0.5) - cosets(i,:))/2, 0);
            aux(i,:) = (aux(i,:)*2 + cosets(i,:));
            dist(i) = norm(x/sqrt(0.5) - aux(i,:))^2;
        end
        index = find(dist==min(dist));
        dec_point = aux(index,:)*sqrt(0.5);
        
        
    case 'E8'   %The Gosset lattice E8. It can be seen as the union of 2 cosets of D8
        cosets = [zeros(1, length(x)); ones(1, length(x))/2];
        dist = zeros(1,2);
        for i=1:size(cosets,1)
            aux(i,:) = lattice_decoding(x - cosets(i,:), 'Dn');
            %aux(i,:) = lattice_decodingc(x - cosets(i,:), 3);
            aux(i,:) = aux(i,:) + cosets(i,:);
            dist(i) = norm(x - aux(i,:))^2;
        end
        index = find(dist==min(dist));
        dec_point = aux(index,:);
end


