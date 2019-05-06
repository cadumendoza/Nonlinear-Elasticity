function [W,S,CC]=NeoHookean_3(C,lambda,mu,icode)

J2=C(1)*C(2)-C(3)*C(3);
J=sqrt(J2);
logJ=log(J);
W = 1/2*lambda*logJ^2 - mu * logJ + 1/2*mu*(C(1)+C(2)-2);

S =[];
CC=zeros(3);

C_inv=[C(2) C(1) -C(3)]/J2;
S= lambda*logJ*C_inv + mu*([1 1 0]-C_inv);
CC(1,1)=lambda*C_inv(1)*C_inv(1) -  ...
    (lambda*logJ-mu)*(C_inv(1)*C_inv(1)+C_inv(1)*C_inv(1));
CC(1,2)=lambda*C_inv(1)*C_inv(2) -  ...
    (lambda*logJ-mu)*(C_inv(3)*C_inv(3)+C_inv(3)*C_inv(3));
CC(1,3)=lambda*C_inv(1)*C_inv(3) -  ...
    (lambda*logJ-mu)*(C_inv(1)*C_inv(3)+C_inv(1)*C_inv(3));
CC(2,2)=lambda*C_inv(2)*C_inv(2) -  ...
    (lambda*logJ-mu)*(C_inv(2)*C_inv(2)+C_inv(2)*C_inv(2));
CC(2,3)=lambda*C_inv(2)*C_inv(3) -  ...
    (lambda*logJ-mu)*(C_inv(2)*C_inv(3)+C_inv(2)*C_inv(3));
CC(3,3)=lambda*C_inv(3)*C_inv(3) -  ...
    (lambda*logJ-mu)*(C_inv(1)*C_inv(2)+C_inv(3)*C_inv(3));
CC(2,1)=CC(1,2);
CC(3,1)=CC(1,3);
CC(3,2)=CC(2,3);
