function [W]=NeoHookean_1(C,lambda,mu,icode)

J2=C(1)*C(2)-C(3)*C(3);
J=sqrt(J2);
logJ=log(J);
W = 1/2*lambda*logJ^2 - mu * logJ + 1/2*mu*(C(1)+C(2)-2);
