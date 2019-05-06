function [Ener,grad_E,Hess_E] = Ener_short(x_short,icode)
global mod1 mesh1 load1 el1

[x_long] = long(x_short);
[Ener,grad_E_l,Hess_E_l] = Energy(x_long,icode);
if icode>1
[grad_E] = short(grad_E_l);
end
if icode == 3
Hess_E = Hess_E_l;
Hess_E(load1.dofCC,:) = [];
Hess_E(:,load1.dofCC) = [];
end