function [Ener,grad_E,Hess_E] = Ener_short(x_short,icode,spring,K,iload)
global mod1 mesh1 load1 el1

[x_long] = long(x_short);

if spring == 1
        x_sp=x_long;
    %     adj=82*ones(41,1);
    %     ADJ=[adj ; -adj];

    %K calculus
    x_ul=(length(load1.dofSpx));
    y_ul=(length(load1.dofSp));
    for i=1:x_ul
        if i==1
            k_sp=K*abs((x_sp(load1.dofSpx(2))-x_sp(load1.dofSpx(1))))/2;
            load1.k_sp(load1.dofSpx(i)+1)=k_sp; 
        elseif i==x_ul
            k_sp=K*abs((x_sp(x_ul)-x_sp(x_ul-2)))/2;
            load1.k_sp(load1.dofSpx(i)+1)=k_sp; 
        else          
            k_sp=K*abs((x_sp(load1.dofSpx(i)+2)-x_sp(load1.dofSpx(i))));
            load1.k_sp(load1.dofSpx(i)+1)=k_sp;           
        end
    end

    %Energy and Grad calculation 
    for i=1:y_ul
        if iload==1
            force_sp=0;
            ener_sp=0;
        else
            force_sp=load1.k_sp(load1.dofSp(i))*abs((mesh1.history_x(iload,load1.dofSp(i))-mesh1.history_x(iload-1,load1.dofSp(i))));
            ener_sp=0.5*load1.k_sp(load1.dofSp(i))*abs((mesh1.history_x(iload,load1.dofSp(i))-mesh1.history_x(iload-1,load1.dofSp(i))))^2;          
        end
        load1.force(load1.dofSp(i))=force_sp;
        load1.enersp(load1.dofSp(i))=ener_sp; 
    end
    
end

[Ener,grad_E_l,Hess_E_l] = Energy(x_long,icode);
for i=1:length(load1.k_sp)
    %Hess_E_l(i,i)=Hess_E_l(i,i)+load1.k_sp(i);
end

if icode>1
[grad_E] = short(grad_E_l);
end
if icode == 3
Hess_E = Hess_E_l;
Hess_E(load1.dofCC,:) = [];
Hess_E(:,load1.dofCC) = [];
load1.k_sp(load1.dofCC)=[];
end

