% This script checks the implementation of the out-of-balance forces 
% (gradient of the potential energy), and the stiffness matrix 
% (the Hessian of the potential energy)

close all
clear all
global mod1 mesh1 load1 el1 undeformed1


example=2;
material=1;
spring=1;   % 1 - with spring, 0 - without
K=0.1;        % Spring constant
[dof_force, dof_disp, lambda, x_eq, CC0, CC1, force, codeLoad]=preprocessing(example,material,spring);
load1.force = force*lambda(1); % include external forces

%Setup the undeformed configuration
precompute;

%define an arbitrary non-homogeneous deformation
%cheching the implementation around the stress-free configuration
%may overlook errors in the geometric stiffness or errors
%compensated by the symmetry of the mesh
x=x_eq;
x(1:2:end) = x(1:2:end) + .5*exp(-2*(x(2:2:end)-.5).^2);
x(2:2:end) = x(2:2:end) + .1*sin(4*x(1:2:end));
%DibujaMalla(mesh1.T,mesh1.x0,x,'r',1)
force_sp=zeros(41,1);
force_sp2=zeros(41,1);
load1.Ensp=0;
load1.Ks=zeros(41,1);
load1.fsp=zeros(328,1);
if spring == 1
    x_sp=x;
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
%compute the energy, its analytical gradient and its analytical Hessian
[Ener,grad_E,Hess_E] = Energy(x,3);

h=1e-8; %for numerical differentiation

%loop over all degrees of freedom
for idof=1:length(x)
   x_=x;
   x_(idof) = x_(idof)+h; % perturb one degree of freedom
   if spring == 1
        x_sp=x_;
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
   [Ener_,grad_E_,Hess_E_] = Energy(x_,3); %compute the perturbed energy and forces
   % 1. Check the gradient
   num_grad=(Ener_-Ener)/h;
   if abs((num_grad-grad_E(idof))/grad_E(idof))>1e-3
       disp('Warning gradient!!')
       disp(idof)
       disp('Analytical')
       disp(grad_E(idof))
       disp('Numerical')
       disp(num_grad)
       disp('relative error')
       disp(abs((num_grad-grad_E(idof))/grad_E(idof)))
   end
  % 2. Check the Hessian
  zumba=Hess_E_(idof,:);
  num_Hess=((grad_E_-grad_E)/h)';
  %num_Hess=zumba;
  if norm(num_Hess-Hess_E(idof,:))/norm(Hess_E(idof,:))>1e-3
      [val index] = max(abs(num_Hess-Hess_E(idof,:)));
      disp('Warning Hessian!!')
      disp(idof)
      disp('Analytical')
      disp(Hess_E(idof,index))
      disp('Numerical')
      disp(num_Hess(index))
      disp('relative error')
      disp(norm(num_Hess-Hess_E(idof,:))/norm(Hess_E(idof,:)))
  end   
end

