function [Ener,grad_E,Hess_E,Ensp] = Ener_short(x_short,icode,spring,K)
global mod1 mesh1 load1 el1

[x_long] = long(x_short);
load1.Ensp=0;
load1.Ks=zeros(41,1);

if spring == 1
    x_sp=x_long;
    force_sp=zeros(41,1);
    load1.dofSpm=load1.dofSp(2:end-1); %degrees of freedom of each y node to assign spring
    force_sp(2:end-1)=-K*0.5*abs((x_sp(load1.dofSpm+1)-x_sp(load1.dofSpm-3))).*(x_sp(load1.dofSpm)-load1.fixedSp(2:end-1));
    force_sp(1)=-K*0.5*abs((x_sp(load1.dofSp(1)+1)-x_sp(load1.dofSp(1)-1))).*(x_sp(load1.dofSp(1))-load1.fixedSp(1));
    force_sp(end)=-K*0.5*abs((x_sp(load1.dofSp(end)-1)-x_sp(load1.dofSp(end)-3))).*(x_sp(load1.dofSp(end))-load1.fixedSp(end));
    load1.Ensp=0.5*force_sp'*(x_sp(load1.dofSp)-load1.fixedSp);
    load1.force(load1.dofSp)=force_sp;
    load1.Ks=[K*0.5*abs((x_sp(load1.dofSp(1)+1)-x_sp(load1.dofSp(1)-1)));K*0.25*abs((x_sp(load1.dofSpm+1)-x_sp(load1.dofSpm-3)));K*0.5*abs((x_sp(load1.dofSp(end)-1)-x_sp(load1.dofSp(end)-3)))];
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