function [W,S,CC]=KsV_3(C,lambda,mu,icode)

% Green-Lagrange strain tensor
E = (C-[1 1 0])/2;
trE = E(1)+E(2);
Esq = E.^2;
trEsq = Esq(1)+Esq(2);
% Energy
W = 0.5*lambda*(trE)^2+mu*trEsq;

% Second Piola-Kirchhoff tensor
S = lambda*trE*[1 1 0]+2*mu*E;

% Constitutive tensor
CC = zeros(3);
CC(1,1) = lambda+2*mu; CC(1,2) = lambda; CC(1,3) = 0;
CC(2,1) = lambda; CC(2,2) = lambda+2*mu; CC(2,3) = 0;
CC(3,1) = 0; CC(3,2) = 0; CC(3,3) = 2*mu;
end