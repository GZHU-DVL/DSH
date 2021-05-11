function [P, theta, sigma2k] = bound_OVE(xk, yk, gamma, P, theta, sigma2k)

%[P, theta, sigma2k] = bound_intersection(xk, yk, gamma, P, theta, sigma2k)
%
%this is the algorithm described in "An optimal volume ellipsoid algorithm
%for parameter set estimation", IEEE Trans. on Automatic Control
%
%Inputs:
%   xk: vector normal to the bounding hyperplanes
%   yk: the current observation multiplied by xk
%   gamma: the independent term in the SME definition
%   P: positive semidefinite matrix which defines the initial ellipsoid
%   theta: center of the initial ellipsoid
%   sigma2k: independent term in the equation of the initial ellipsoid
%
%Outputs:
%   P: positive semidefinite matrix of the updated ellipsoid
%   theta: center of the updated ellipsoid
%   sigma2k: independent term in the equation of the updated ellipsoid


%P
dimensions = size(P,1);
deltak = yk - xk'*theta;
Gk = xk'*P*xk;

alphak = (deltak + gamma)/sqrt(Gk);
Betak = gamma/sqrt(Gk);
% alphak
% [-1, 1+2*Betak]
%[alphak-2*Betak, alphak]
if (alphak<-1) || ((alphak-2*Betak)>1)
    %disp('***** No cut *****') %no intersection
    P = zeros(dimensions);
    theta = zeros(dimensions,1);
    sigma2k = 1;
    return
end
if (alphak>1) && ((alphak-2*Betak)<-1)
    %disp('No update')  %the hyperplanes do not cut the current ellipsoid
    return    %no update
end
if alphak > 1
    Betak = Betak -(alphak-1)/2;
    alphak = 1;
end
if (2*Betak-alphak) > 1
    Betak = (1 + alphak)/2;
end
%[alphak-2*Betak alphak]
if alphak==Betak
    tauk = 0;
    lambdak = dimensions*(1-Betak^2)/(dimensions-1);
    sigmak = (1-dimensions*Betak^2)/(1-Betak^2);
else
    aux_tau = linspace(alphak-2*Betak + Betak/10000, alphak - Betak/10000, 10000);
    %             [aux_tau(1), aux_tau(end)]
    poly_vector = [(dimensions+1), ((1+alphak)*(alphak-2*Betak+1)/(Betak-alphak) + ...
        2*(dimensions*(Betak-alphak)+1)), dimensions*alphak*(alphak-2*Betak) + 1];
    roots_vector = roots(poly_vector);
    if (roots_vector(1)>aux_tau(1))&(roots_vector(1)<aux_tau(end))
        tauk = roots_vector(1);
    else
        tauk = roots_vector(2);
    end
    %             equation = (dimensions+1).*aux_tau.^2 + ((1+alphak)*(alphak-2*Betak+1)/(Betak-alphak) + 2*(dimensions*(Betak-alphak)+1)).*aux_tau ...
    %                 + dimensions*alphak*(alphak-2*Betak) + 1;
    %             indexmin = find(abs(equation)==min(abs(equation)), 1);
    %             figure(10), plot(aux_tau, equation), grid on, pause
    lambdak = ((tauk + 1)^2*(Betak - alphak) - tauk*(1 + alphak)*(2*Betak - alphak - 1))/(tauk + Betak - alphak);
    sigmak = -tauk/(Betak - alphak);
end
if sigmak < 0
    %disp('No update')  %it is not possible to obtain another bounding ellipsoid with smaller volume than the current one
    return
end
theta = theta + (tauk*P*xk)/sqrt(Gk);
P = lambdak*(P - sigmak*(P*xk*xk'*P)/Gk);
