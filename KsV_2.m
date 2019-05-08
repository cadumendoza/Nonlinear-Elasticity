function [W,S]=KsV_2(C,lambda,mu,icode)

% Green-Lagrange strain tensor
E = (C-[1 1 0])/2;
trE = E(1)+E(2);
Esq = E.^2;
trEsq = Esq(1)+Esq(2);

% Energy
W = 0.5*lambda*(trE)^2+mu*trEsq;

% Second Piola-Kirchhoff tensor
S = lambda*trE*[1 1 0]+2*mu*E;
end