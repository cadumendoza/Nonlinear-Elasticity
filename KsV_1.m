function [W]=KsV_1(C,lambda,mu,icode)

%tr(E)
trE=0.5*(C(1)+C(2))-1;
%tr(E^2)
E11=0.5*(C(1)-1);
E12=0.5*(C(3));
E22=0.5*(C(2)-1);

trE2=E11*E11+2*E12*E12+E22*E22;

W = 1/2*lambda*trE^2+mu*trE2;
