function [W,S,CC]=transv_isotr_3(C,c0,c1,kappa,mu,N)

J2 = C(1)*C(2)-C(3)*C(3); % Jacobian square
J = sqrt(J2);
lnJ = log(J);

trC = C(1)+C(2); % trace of C
invC = [C(2) C(1) -C(3)]/J2; % Inverse of C
Iv = [1 1 0]; % Voigt notation of identity matrix
N_fib = [N(1)*N(1), N(2)*N(2), N(1)*N(2)];

% Derivatives
dC11=-1/2*(invC(1)*invC(1)+invC(1)*invC(1));
dC22=-1/2*(invC(2)*invC(2)+invC(2)*invC(2));
dC33=-1/2*(invC(1)*invC(2)+invC(3)*invC(3));
dC12=-1/2*(invC(3)*invC(3)+invC(3)*invC(3));
dC13=-1/2*(invC(1)*invC(3)+invC(1)*invC(3));
dC23=-1/2*(invC(2)*invC(3)+invC(2)*invC(3));

% Volumetric response of the material
G = 0.25*(J2-1-2*lnJ);
% Fourth invariant
I4 = C(1)*N(1)*N(1)+C(2)*N(2)*N(2)+2*C(3)*N(1)*N(2);
% Energy
W = mu*0.5*(trC-2)-mu*lnJ+kappa*G+c0*(exp(c1*(sqrt(I4)-1)^4)-1);
% Second Piola-Kirchhoff tensor
S = mu*(Iv-invC)+kappa/2*(J2*invC-invC)+4*c0*c1*...
((sqrt(I4)-1)^3/sqrt(I4)*exp(c1*(sqrt(I4)-1)^4))*N_fib;

% Constitutive tensor
value=(8*c0*c1*((sqrt(I4)-1)^2)/I4)*exp(c1*((sqrt(I4)-1)^4))*...
(2*c1*((sqrt(I4)-1)^4)+(3/2)-((sqrt(I4)-1)/(2*sqrt(I4))));

CC(1,1)=(-2*mu + kappa*(J2-1))*dC11 + kappa*J2*invC(1)*invC(1) +...
value*N_fib(1)*N_fib(1);

CC(2,2)=(-2*mu + kappa*(J2-1))*dC22 + kappa*J2*invC(2)*invC(2) +...
 value*N_fib(2)*N_fib(2);

CC(3,3)=(-2*mu + kappa*(J2-1))*dC33 + kappa*J2*invC(3)*invC(3) +...
value*N_fib(3)*N_fib(3);

CC(1,2)=(-2*mu + kappa*(J2-1))*dC12 + kappa*J2*invC(2)*invC(1) +...
 value*N_fib(1)*N_fib(2);

CC(1,3)=(-2*mu + kappa*(J2-1))*dC13 + kappa*J2*invC(1)*invC(3) +...
    value*N_fib(1)*N_fib(3);

CC(2,3)=(-2*mu + kappa*(J2-1))*dC23 + kappa*J2*invC(2)*invC(3) +...
value*N_fib(2)*N_fib(3);

CC(2,1)=CC(1,2);
CC(3,1)=CC(1,3);
CC(3,2)=CC(2,3);
end