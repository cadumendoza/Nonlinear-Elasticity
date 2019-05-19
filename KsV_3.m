function [W,S,CC]=KsV_3(C,lambda,mu,icode)

%---------------------------------------%
% Kirchhoff-Saint Venant Material Model %
%---------------------------------------%

% Trace of Green-Lagrange Strain Tensor E
% tr(E)
trE=0.5*(C(1)+C(2))-1;

% Trace of the square of Green-Lagrange Strain Tensor E
% tr(E^2)
E11=0.5*(C(1)-1);
E12=0.5*(C(3));
E22=0.5*(C(2)-1);

trE2=E11*E11+2*E12*E12+E22*E22;

% Energy Function
W = 1/2*lambda*trE^2+mu*trE2;

% Second Piola-Kirchhoff Stress Tensor
S =[];
S= lambda*trE*[1 1 0] + 2*mu*[E11 E22 E12];

% Constitutive Elastic Tensor
CC=zeros(3);

CC(1,1)=lambda+2*mu;
CC(1,2)=lambda;
CC(2,1)=lambda;
CC(2,2)=lambda+2*mu;
CC(3,3)=mu;