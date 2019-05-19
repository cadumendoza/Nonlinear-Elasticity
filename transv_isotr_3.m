function [W,S,CC]=transv_isotr_3(C,c0,c1,kappa,mu,N)

%---------------------------------------%
% Transversely Isotropic Material Model %
%---------------------------------------%

J2 = C(1)*C(2)-C(3)*C(3); % Jacobian square
J = sqrt(J2);
lnJ = log(J);

% Trace of Right Cauchy-Green Strain Tensor C
% tr(C)
trC = C(1)+C(2); % trace of C
invC = [C(2) C(1) -C(3)]/J2; % Inverse of C
Iv = [1 1 0]; % Voigt notation of identity matrix
%Fibers direction
N_f = [N(1)*N(1), N(2)*N(2), N(1)*N(2)];

% Derivatives
DC11=-1/2*(invC(1)*invC(1)+invC(1)*invC(1));
DC22=-1/2*(invC(2)*invC(2)+invC(2)*invC(2));
DC33=-1/2*(invC(1)*invC(2)+invC(3)*invC(3));
DC12=-1/2*(invC(3)*invC(3)+invC(3)*invC(3));
DC13=-1/2*(invC(1)*invC(3)+invC(1)*invC(3));
DC23=-1/2*(invC(2)*invC(3)+invC(2)*invC(3));

% Fourth invariant
I4 = C(1)*N(1)*N(1)+C(2)*N(2)*N(2)+2*C(3)*N(1)*N(2);

% Volumetric response of the material
Gfunc = 0.25*(J2-1-2*lnJ);

% Energy Function
W = mu*0.5*(trC-2)-mu*lnJ+kappa*Gfunc+c0*(exp(c1*(sqrt(I4)-1)^4)-1);

% Second Piola-Kirchhoff Stress Tensor
S = mu*(Iv-invC)+kappa/2*(J2*invC-invC)+4*c0*c1*((sqrt(I4)-1)^3/sqrt(I4)*exp(c1*(sqrt(I4)-1)^4))*N_f;

% Constitutive Elastic Tensor
value=(8*c0*c1*((sqrt(I4)-1)^2)/I4)*exp(c1*((sqrt(I4)-1)^4))*(2*c1*((sqrt(I4)-1)^4)+(3/2)-((sqrt(I4)-1)/(2*sqrt(I4))));

CC(1,1)=(-2*mu + kappa*(J2-1))*DC11 + kappa*J2*invC(1)*invC(1) + value*N_f(1)*N_f(1);

CC(2,2)=(-2*mu + kappa*(J2-1))*DC22 + kappa*J2*invC(2)*invC(2) + value*N_f(2)*N_f(2);

CC(3,3)=(-2*mu + kappa*(J2-1))*DC33 + kappa*J2*invC(3)*invC(3) + value*N_f(3)*N_f(3);

CC(1,2)=(-2*mu + kappa*(J2-1))*DC12 + kappa*J2*invC(2)*invC(1) + value*N_f(1)*N_f(2);

CC(1,3)=(-2*mu + kappa*(J2-1))*DC13 + kappa*J2*invC(1)*invC(3) + value*N_f(1)*N_f(3);

CC(2,3)=(-2*mu + kappa*(J2-1))*DC23 + kappa*J2*invC(2)*invC(3) + value*N_f(2)*N_f(3);

CC(2,1)=CC(1,2);
CC(3,1)=CC(1,3);
CC(3,2)=CC(2,3);
end