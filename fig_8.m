%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dither estimation using averaging of the modulo-lambda reduced observations (KMA case)
%It is the optimal estimator in infinite dimensions when the embedding lattice is the 
% optimal lattice quantizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
%parameters
lattice = 'E8'   %the type of lattice considered in the experiment
%dimensions = 2;  %number of dimensions of the considered lattice
num_iter = 5;   %number of realizations
alpha = [0.3, 0.5, 0.8];  %distortion compensation parameter
No = [10:10:100];  %number of observations
dimensions = 8;
vol = 1;
GE8 = 0.071682;
Delta = sqrt(vol^(-2/dimensions)/(GE8*12));
%No = [2:2:10,15,20,25,30];
%No = [10,20:10:100];
debug = 0;   %indicates whether debugging information is plotted or not
rnstate = 0;    %initial state of the random number generator
%the scaling factor Delta and the covering radius are computed
% switch lattice
%     case 'cubic'
%         Delta = 1
%         lattice_type = 0;
%         cov_radius = ones(1,dimensions)*Delta*sqrt(dimensions)/2  %cubic lattice
%     case 'hexagonal'
%         dimensions = 2;
%         Ghexagonal = 0.080188
%         M = [1 0; 1/2 sqrt(3)*1/2];
%         vol = sqrt(det(M*M'));
%         Delta = sqrt((1/12)/(Ghexagonal*vol))
%         lattice_type = 1;
%         %cov_radius = Delta/sqrt(3); %hexagonal lattice
%         cov_radius = [Delta/sqrt(3), Delta/2] %hexagonal lattice
%     case 'sphere'
%         Gsphere = (gamma(0.5*dimensions+1)^(2/dimensions)/((dimensions+2)*pi))
%         volsphere = (pi^(dimensions/2)/gamma(dimensions/2 + 1))*(1/2)^(dimensions);  %volume of the sphere of radius 0.5;
%         Delta = sqrt(volsphere^(-2/dimensions)/(Gsphere*12));
%         cov_radius = ones(1,dimensions)*Delta/2
%         lattice_type = 2;
%     case 'Dn'
%         GDn = (1/12)*(dimensions==2) + 0.0078543*(dimensions==3) + 0.076603*(dimensions==4) + 0.075786*(dimensions==5)
%         vol = 2;
%         Delta = sqrt(vol^(-2/dimensions)/(GDn*12));
%         cov_radius = ones(1,dimensions)*(Delta*(dimensions==2) + Delta*1*(dimensions==3) + ...
%             (Delta*sqrt(dimensions/2)/sqrt(2))*(dimensions>3))
%         lattice_type = 3;
%     case 'Dndual'
%         GDndual = 0.078543*(dimensions==3) + 0.076603*(dimensions==4) + 0.075625*(dimensions==5)
%         vol = 0.5;
%         Delta = sqrt(vol^(-2/dimensions)/(GDndual*12))
%         cov_radius = ones(1,dimensions)*((Delta*sqrt(5/3)*sqrt(3)/4)*(dimensions==3) + ...
%             (Delta*sqrt(dimensions)/(2*sqrt(2)))*((dimensions>=4)&&(mod(dimensions,2)==0)) + ...
%             (Delta*sqrt(2*dimensions-1)/4)*((dimensions>=5)&&(mod(dimensions,2)~=0)))
%         lattice_type = 4;
%     case 'E7dual'
%         dimensions = 7;
%         GE7dual = 0.073116
%         vol = sqrt(0.5);
%         Delta = sqrt(vol^(-2/dimensions)/(GE7dual*12));
%         cov_radius = ones(1, dimensions)*sqrt(7/8)*Delta
%         lattice_type = 6;
%     case 'E8'
%         dimensions = 8;
%         GE8 = 0.071682;
%         vol = 1;
%         Delta = sqrt(vol^(-2/dimensions)/(GE8*12))
%         cov_radius = ones(1,dimensions)*Delta*1;
%         lattice_type = 5;
% end
randn('state', 0)
error=zeros(10,3);
for k=1:3
for i=1:num_iter
    %observations uniformly distributed around the origin
    obs_iter = (1-alpha(k)).*rand_obs(No(end), Delta, lattice, dimensions)';
    for j=1:length(No)           
        estimate = mean(obs_iter(:,1:j), 2);
        serror(j,i) = (norm(estimate - zeros(dimensions,1))^2)/dimensions;
    end
end
error(:,k)=mean(serror,2);
end
semilogy(No, error(:,2)','-.k', No, error(:,3)','-^r', 'Linewidth', 1);
legend('SME-Based','Proposed');
xlabel('Number of observations','FontSize',18,'FontName','Times New Roman');
ylabel('MSE','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',18,'FontName','Times New Roman');


