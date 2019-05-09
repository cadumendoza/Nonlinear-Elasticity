function [Ener,grad_E,Hess_E] = Ener_short(x_short,icode,spring,K)
global mod1 mesh1 load1 el1

[x_long] = long(x_short);

if spring == 1
    x_sp=x_long;
    force_sp=-K*abs((x_sp(load1.dofSp+2)-x_sp(load1.dofSp))).*x_sp(load1.dofSp)
    load1.force(load1.dofSp)=force_sp;
end
[Ener,grad_E_l,Hess_E_l] = Energy(x_long,icode);
if icode>1
[grad_E] = short(grad_E_l);
end
if icode == 3
Hess_E = Hess_E_l;
Hess_E(load1.dofCC,:) = [];
Hess_E(:,load1.dofCC) = [];
end