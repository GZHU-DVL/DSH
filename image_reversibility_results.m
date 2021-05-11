%reversibility of a real host signal (host estimation) based on the secret dither estimate
clear all, close all, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lattice = 'hexagonal';
alpha=0.75;
image_file = 'couple.bmp' ; %we select an image to be marked (leave empty for Gaussian observations)
%image_file = []  %leave empty for Gaussian observations
indices = [2 3 4 9 10 11 12 17 18 19 25 26];    %indices (in a 8x8 block) of the DCT coefficients to be watermarked (in case an image will be used)
DWR = 30;
switch lattice
    case 'hexagonal'
        method = 'dcdm_generic';
        n = 2;
        g = [1 2]';
        sigmax=10;
        M = [0, 1; sqrt(3)/2, 1/2];
        vol = sqrt(det(M*M'));
        Ghexagonal = 0.080188;  %normalized second order moment
        Delta = sqrt((sigmax^2*10^(-DWR/10)/(alpha^2))/(Ghexagonal*vol));
end
    im = im2double(imread(image_file));
    [a,b,c] = size(im);
    if a ~= b
        disp('error: The input image must be square');
        return
    end;
    if size(im, 3)~=1   %if the image is not in grayscale, it is converted
        im = rgb2gray(im);
    end
    figure(1), imshow(im); % title('original image');        
    d = 8;  %size of the DCT blocks 
    nblocks = floor(size(im, 2)/d); %number of DCT blocks    
    x = zeros(n, nblocks^2);
    nobs = 1;
    tx = zeros(d^2, nblocks^2);
    for i=0:nblocks-1   %vertical index
        for q=0:nblocks-1   %horizontal index
            temp = im((i*d)+1 : ((i+1)*d), (q*d)+1 : ((q+1)*d));
            z = dct2(temp);
            z = reshape(z, d^2, 1);
            tx(:,nobs) = z; %64*16384
            x(:,nobs) = z(indices(1:n),1);%2*16384
            nobs = nobs + 1;
        end
    end
    normfactor = sqrt(40);
    x = normfactor*x;
%first, the whole set of observations (host image) is watermarked with the true dither signal
p=9;
switch  method
    case 'dcdm_generic'
        [obs2, host, cosets, alphabet, message, watermark, Delta, alpha, cov_radius, dither, DWR] = ...
            dcdm_generic(x, Delta, alpha, lattice, n, p, g, 0, DWR);
end
%reconstruction of the watermarked image in the spatial domain
nobs = 1;
im_wat = im;
for i=0:nblocks-1   %vertical index
    for q=0:nblocks-1   %horizontal index
        tx(indices(1:n),nobs) = obs2(:,nobs)/normfactor;
        im_wat((i*d)+1 : ((i+1)*d), (q*d)+1 : ((q+1)*d)) = idct2(reshape(tx(:,nobs), d, d));
        nobs = nobs + 1;
    end
end
figure(2), imshow(im_wat); %title('watermarked image');
%psnr=10*log10(size(im,1)*size(im,2)*255.^2) - 10*log10(sum(sum((im2uint8(abs(im - im_wat)).^2))))    %PSNR
%now, reversibility is performed with the most likely estimated dither vector
dither_estimate=dither;
dec_message = dcdm_decoding(lattice, Delta, dither_estimate, message, cosets, obs2);
est_host = reverse_dcdm(lattice, Delta, alpha, dither_estimate, dec_message', cosets, obs2);
    %reconstruction of the estimated host in the spatial domain
    nobs = 1;
    im_est = im;
    for i=0:nblocks-1   %vertical index
        for q=0:nblocks-1   %horizontal index
            tx(indices(1:n),nobs) = est_host(:,nobs)/normfactor;
            im_est((i*d)+1 : ((i+1)*d), (q*d)+1 : ((q+1)*d)) = idct2(reshape(tx(:,nobs), d, d));
            nobs = nobs + 1;
        end
    end
    fprintf('PSNR: %f\n', 10*log10(a*b*255^2) - 10*log10(sum(sum((255*(im - im_est)).^2))))    %PSNR    
    figure(3), imshow(im_est); %title('Estimated host')
    figure(4), stem(message(1:100) - dec_message(1:100));

