function [W,S]=transv_isotr_2(C,c0,c1,kappa,mu,N)

J2 = C(1)*C(2)-C(3)*C(3); % Jacobian square
J = sqrt(J2);
lnJ = log(J);

trC = C(1)+C(2); % trace of C
invC = [C(2) C(1) -C(3)]/J2; % Inverse of C
Iv = [1 1 0]; % Voigt notation of identity matrix
N_fib = [N(1)*N(1), N(2)*N(2), N(1)*N(2)];

% Volumetric response of the material
G = 0.25*(J2-1-2*lnJ);
% Fourth invariant
I4 = C(1)*N(1)*N(1)+C(2)*N(2)*N(2)+2*C(3)*N(1)*N(2);
% Energy
W = mu*0.5*(trC-2)-mu*lnJ+kappa*G+c0*(exp(c1*(sqrt(I4)-1)^4)-1);
% Second Piola-Kirchhoff tensor
S = mu*(Iv-invC)+kappa/2*(J2*invC-invC)+4*c0*c1*...
((sqrt(I4)-1)^3/sqrt(I4)*exp(c1*(sqrt(I4)-1)^4))*N_fib;

end