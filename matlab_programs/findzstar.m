%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: findzstar.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: reports the difference in profits between the modern and
% traditional sectors for a firm with productivity zstar in the analytic
% version of the model from Section 3 of the paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function diff = findzstar(zstar, lambda, eta, gamma, v, Am, vstar, Omega, R)

N =1; 
Zbart = lambda/(lambda-1)*(1-zstar^(1-lambda))/(1-zstar^(-lambda)); 
Zbarm = (lambda/(lambda-1))*zstar; 
Nm_N = zstar^(-lambda); %Fraction of modern firms 
Nm = Nm_N*N; 
Nt = (1-Nm_N)*N; 

Kt = Nt*Zbart*(eta/R)^(1/(1-eta)); 
Km0 = (1-gamma)*Nm*Zbarm*(v*eta*Am/R).^(1/(1-eta)); 

Km1_l = gamma*Nm*Zbarm*(eta*Am./(R + R  + 1-v)).^(1/(1-eta)); 
Km1_h = gamma*Nm*Zbarm*(v*eta*Am/R).^(1/(1-eta)); 
Ks1_l = Km1_l;
Ks1_h = gamma*Nm*Zbarm*((1-v)*eta*Am./(R + 1-v)).^(1/(1-eta)); 

if (v < vstar)
    Km1 = Km1_l; 
    Ks1 = Ks1_l; 
else
    Km1 = Km1_h; 
    Ks1 = Ks1_h; 
end

yhatt = (Kt/(Nt*Zbart))^eta; 
yhatm1 = Am^(1-eta)*(v*(Km1/(gamma*Nm*Zbarm))^eta + (1-v)*(Ks1/(gamma*Nm*Zbarm))^eta);
yhatm0 = Am^(1-eta)*v*(Km0/((1-gamma)*Nm*Zbarm))^eta; 

diff = zstar*(1-eta)*(gamma*yhatm1 + (1-gamma)*yhatm0  - yhatt) - Omega; 
end