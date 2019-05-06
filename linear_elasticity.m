function [u,Reaction,delta] = linear_elasticity(codeLoad,fact,example,CC1,dof_force,dof_disp)
global mod1 mesh1 load1 el1 undeformed1


x_eq=mesh1.x0;
x=x_eq;
switch example
   case {1, 2}
      x(1:2:end)=x_eq(1:2:end)*fact;
end

x_BC = x(load1.dofCC);
x0_BC=x_eq(load1.dofCC);
u_BC=x_BC-x0_BC;

[E_eq,grad_eq, K0] = Energy(x_eq,3);
L=zeros(size(K0,1),load1.ndofCC);
Zero_M=zeros(load1.ndofCC,load1.ndofCC);
for i=1:load1.ndofCC
    L(load1.dofCC(i),i) = 1;
end
K_L = [K0,L;L',Zero_M];
f=[zeros(size(K0,1),1);u_BC] + [load1.force; zeros(size(u_BC))];

sol=K_L\f;
u=sol(1:size(K0,1));
Reac = sol(size(K0,1)+1:end);
nn=length(Reac);


switch example
    case {1, 2}
    Reaction = sum(Reac(1:nn/4));
    delta = u(2*CC1(1)-1);
    case {0, 3}
%    Reaction = sum(Reac(1:nn/3));
    Reaction =sum(load1.force); 
    delta = u(2*CC1(1)-1);
    case {4, 5}
        Reaction =mean(load1.force(dof_force)); 
        delta = mean(u(dof_disp));
    otherwise
        disp('Case not implemented')
end

