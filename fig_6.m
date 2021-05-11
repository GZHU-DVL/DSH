%volume in a lattice through Monte Carlo simulations.
%The intersection region is computed by Montecarlo. The empirical volume is
% the average of log(volume) over num_iter realizations
clear all; close all;
clc
lattice = 'hexagonal';   
dimension = 2;  %number of dimensions of the considered lattice
N=5000;  % the search points in voronoi region
samples = N;    %number of samples uniformly generated on the bounding region
num_iter = 1;   %number of realizations
alph =[0.5:0.02:1.0];  %distortion compensation parameter
debug = 0;   %indicates whether debugging information is plotted or not (works only for 2 dimensions)
rnstate = 0;    %initial state of the random number generator
p=9;
alpha=0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first, the scaling factor Delta and the covering radius are computed
Ghexagonal = 0.080188;
M = [1 0; 1/2 sqrt(3)*1/2];
vol = sqrt(det(M*M'));
Delta = sqrt((1/12)/(Ghexagonal*vol));
lattice_type = 1;   
Ghex = [0, Delta; Delta*sqrt(3)/2, Delta/2];
cov_radius = [Delta/sqrt(3), Delta/2]; %hexagonal lattice
g=[2 3]';
dither=[0 0]';
cosets = construction_a(p, g, Ghex, 0);
rand('state', rnstate);  %the seed of the random generator is fixed
%this is the MonteCarlo algorithm
volume_int = zeros(length(alpha),num_iter);
 theta=0;
 R=[cos(theta), -sin(theta); sin(theta), cos(theta)];
 num_obs=30; % Monte Carlo number to approximate SER (symbol error rate)
 No=num_obs;
 a=[(1-p.^(-1/dimension)):0.02:1];
 s_v=sqrt(3)*Delta/2;
% for i=1:num_iter
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%% uniform search beginning
% obs =R*rand_obs(num_obs, Delta, lattice, dimension)';
% message = floor(rand(1, No)*p); %sequence of symbols (p-ary, equiprobable) to be embedded
% alphabet = 0:p-1;
% message = mod(message - message(1), length(alphabet)); 
% % watermark embedding
% for index_obs = 1:No
%     host = obs(:,index_obs);
%     coset_index = message(index_obs)+1;    
%     quant = lattice_decoding((host' - cosets(:,coset_index)' - dither')/Delta , lattice_type)*Delta + ...
%         cosets(:,coset_index)' + dither'; % dither(:,index_obs)' take place of dither by wyg
%     watermark(:, index_obs) = alpha*(quant' - host);
%     obs(:,index_obs) = host + watermark(:, index_obs);
% end
% point=N;
% vol_voronoi=(sqrt(3)*Delta.^2)/2;
% vol_point=vol_voronoi/point;
% step_search=sqrt(vol_point);
% % four linear equations listed as follows from the 1- to 4-th quadrant
% % respectively, y=-sqrt(3)*x+Delta; y=sqrt(3)*x+Delta; y=-sqrt(3)*x-Delta;
% % y=sqrt(3)*x-Delta
% pn=1;
% est_dither=zeros(2,point);
% %search divided into three regions
% for x=-Delta/sqrt(3):step_search:(-Delta/(2*sqrt(3))-step_search)
%     for y=(-sqrt(3)*x-Delta):step_search:(sqrt(3)*x+Delta)
%         est_dither(:,pn)=[x, y]';
%         % decoding using the estimated dither
%         dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
%         % compute the symbol error rate (SER)
%        ser(pn)=sum(dec_message~=message)/No;
%         pn=pn+1;
%     end
% end
% for x=(-Delta/(2*sqrt(3))):step_search:(Delta/(2*sqrt(3)))
%     for y=(-Delta/2):step_search:(Delta/2)
%         est_dither(:,pn)=[x, y]';
%         % decoding using the estimated dither
%         dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
%         % compute the symbol error rate (SER)
%        ser(pn)=sum(dec_message~=message)/No;
%         pn=pn+1;
%     end
% end
% for x=(Delta/(2*sqrt(3))+step_search):step_search:(Delta/sqrt(3))
%     for y=(sqrt(3)*x-Delta):step_search:(-sqrt(3)*x+Delta)
%         est_dither(:,pn)=[x, y]';
%         % decoding using the estimated dither
%         dec_message = dcdm_decoding(lattice, Delta, est_dither(:,pn), message, cosets, obs);
%         % compute the symbol error rate (SER)
%        ser(pn)=sum(dec_message~=message)/No;
%         pn=pn+1;
%     end
% end
% index_u=find(ser==0);
% U_dither=est_dither(:,index_u); 
% %%%% uniform search end
% 
% %%%% randon search beginning
% est_dither_r=rand_obs(N, Delta, lattice, dimension)';
% for pn=1:N
%     dec_message = dcdm_decoding(lattice, Delta, est_dither_r(:,pn), message, cosets, obs);
%         % compute the symbol error rate (SER)
%     ser(pn)=sum(dec_message~=message)/No;
% end
% index_r=find(ser==0);
% R_dither=est_dither_r(:,index_r); 
% %%%% randon search end
% % 
%     true_up_limits_u = min(1,U_dither(1:length(index_u)));
%     true_low_limits_u = max(2, U_dither(1:length(index_u)));
%    
%     true_up_limits_r = min(1, R_dither(1:length(index_u)));
%     true_low_limits_r = max(2, R_dither(1:length(index_u)));
%     %Now, the volume of the set of tolerable estimates is calculated
%     bounding_volume_u(1,i) = prod(true_up_limits_u-true_low_limits_u); %this is the volume of the region that bounds the intersection
%     bounding_volume_r(1,i)  = prod(true_up_limits_r-true_low_limits_r );  %and this is the volume of the intersection
% end
figure;
vol_t=(a-ones(1,size(a,2))+repmat(p.^(-1/dimension), 1, size(a,2)))*s_v;
plot(a, vol_t, 'r-.','Linewidth', 3);
p=4;
a=[(1-p.^(-1/dimension)):0.02:1];
vol_r=(a-ones(1,size(a,2))+repmat(p.^(-1/dimension), 1, size(a,2)))*s_v;
hold on;
plot(a, vol_t, 'r-.','Linewidth', 3);
plot(a, vol_r,'b-x' ,'Linewidth', 1);
legend('Theoretical','Estimated');
xlabel('$\alpha$','Interpreter','LaTex','FontSize',18,'FontName','Times New Roman');
ylabel('$\textrm{vol}[\mathcal{S}(\hat{\textbf{t}})]$','Interpreter','LaTex','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',18,'FontName','Times New Roman');
hold off;



