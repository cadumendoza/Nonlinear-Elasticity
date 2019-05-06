% This script checks the implementation of the out-of-balance forces 
% (gradient of the potential energy), and the stiffness matrix 
% (the Hessian of the potential energy)

close all
clear all
global mod1 mesh1 load1 el1 undeformed1


example=0;
material=3;
[dof_force, dof_disp, lambda, x_eq, CC0, CC1, force, codeLoad]=preprocessing(example,material);
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
DibujaMalla(mesh1.T,mesh1.x0,x,'r',1)

%compute the energy, its analytical gradient and its analytical Hessian
[Ener,grad_E,Hess_E] = Energy(x,3);

h=1e-8; %for numerical differentiation

%loop over all degrees of freedom
for idof=1:length(x)
   x_=x;
   x_(idof) = x_(idof)+h; % perturb one degree of freedom
   [Ener_,grad_E_] = Energy(x_,2); %compute the perturbed energy and forces
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
  num_Hess=(grad_E_-grad_E)/h;
  num_Hess=num_Hess';
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

