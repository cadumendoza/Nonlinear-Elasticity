function [Ener,grad_E,Hess_E,Ensp] = Ener_short(x_short,icode,spring,K)
global mod1 mesh1 load1 el1

[x_long] = long(x_short);
force_sp=zeros(41,1);
force_sp2=zeros(41,1);
if spring == 1
    x_sp=x_long;
    for i=1:length(load1.dofSp)
        if x_sp(load1.dofSp(i))>0
            x_sp(load1.dofSp(i))=0;
        end
    end
    for i=1:length(load1.dofSp2)
        if x_sp(load1.dofSp2(i))<1
            x_sp(load1.dofSp2(i))=0;
        end
    end
    load1.dofSpm=load1.dofSp(2:end-1);
    force_sp(2:end-1)=-K*0.5*abs((x_sp(load1.dofSpm+1)-x_sp(load1.dofSpm-3))).*(x_sp(load1.dofSpm)-load1.fixedSp(2:end-1));
    force_sp(1)=-K*0.5*abs((x_sp(load1.dofSp(1)+1)-x_sp(load1.dofSp(1)-1))).*(x_sp(load1.dofSp(1))-load1.fixedSp(1));
    force_sp(end)=-K*0.5*abs((x_sp(load1.dofSp(end)-1)-x_sp(load1.dofSp(end)-3))).*(x_sp(load1.dofSp(end))-load1.fixedSp(end));
    load1.dofSpm2=load1.dofSp2(2:end-1);
    force_sp2(2:end-1)=-K*0.5*abs((x_sp(load1.dofSpm2+1)-x_sp(load1.dofSpm2-3))).*(x_sp(load1.dofSpm2)-load1.fixedSp2(2:end-1));
    force_sp2(1)=-K*0.5*abs((x_sp(load1.dofSp2(1)+1)-x_sp(load1.dofSp2(1)-1))).*(x_sp(load1.dofSp2(1))-load1.fixedSp2(1));
    force_sp2(end)=-K*0.5*abs((x_sp(load1.dofSp2(end)-1)-x_sp(load1.dofSp2(end)-3))).*(x_sp(load1.dofSp2(end))-load1.fixedSp2(end));
    load1.Ensp=0.5*force_sp'*(x_sp(load1.dofSp)-load1.fixedSp)+0.5*force_sp2'*(x_sp(load1.dofSp2)-load1.fixedSp2);
    load1.fsp(load1.dofSp)=-force_sp;
    load1.fsp(load1.dofSp2)=-force_sp2;
    load1.Ks=[K*0.5*abs((x_sp(load1.dofSp(1)+1)-x_sp(load1.dofSp(1)-1)));K*0.5*abs((x_sp(load1.dofSpm+1)-x_sp(load1.dofSpm-3)));K*0.5*abs((x_sp(load1.dofSp(end)-1)-x_sp(load1.dofSp(end)-3)))];
    load1.Ks2=[K*0.5*abs((x_sp(load1.dofSp2(1)+1)-x_sp(load1.dofSp2(1)-1)));K*0.5*abs((x_sp(load1.dofSpm2+1)-x_sp(load1.dofSpm2-3)));K*0.5*abs((x_sp(load1.dofSp2(end)-1)-x_sp(load1.dofSp2(end)-3)))];
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