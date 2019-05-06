function [W]=transv_isotr_1(C,c0,c1,kappa,mu,N)

J2 = C(1)*C(2)-C(3)*C(3); % square of the Jacobian
J = sqrt(J2);
lnJ = log(J);
trC = C(1)+C(2); % trace of C

% Volumetric response of the material
G = 0.25*(J2-1-2*lnJ);
% Fourth invariant
I4 = C(1)*N(1)*N(1)+C(2)*N(2)*N(2)+2*C(3)*N(1)*N(2);
% Energy
W = mu*0.5*(trC-2)-mu*lnJ+kappa*G+c0*(exp(c1*(sqrt(I4)-1)^4)-1);
end