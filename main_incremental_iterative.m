close all
clear all
global mod1 mesh1 load1 el1 undeformed1

% 0: upsetting of a block, dead load
% 1: upsetting of a block, imposed displacements
% 2: compression of a slender beam, imposed displacements
% 3: compression of a slender beam, dead load
% 4: arch, dead load at center of the arch
% 5: arch, dead load near the supports
example=2;
material=1;
spring=1;   % 1 - with spring, 0 - without
K=0.01;        % Spring constant
[dof_force, dof_disp, lambda, x_eq, CC0, CC1, force, codeLoad]=preprocessing(example,material,spring);

%Equilibrate
options.n_iter_max=80;
options.tol_x=1.e-6;
options.tol_f=1.e-6;
options.info=3;
options.method=1; %0: vanilla Newton-Rapshon, 1: Newton-Rapshon
%11: Modified NR, 2: L-BFGS, 3: Conjugate Gradient
options.linesearch=1; % 0: off, 1: on. For method 3, automatically on.
options.CGv = 2; %Method 3 is needed  (1: Fletcher-Reeves, 2: Polak-Ribiere)
options.L=8; %Method 2 is needed: value of L for L-BFGS
options.prec=1; % Preconditioner for methods 2&3 (0: No Active, 1: Active)


% Options for Line Search
options.n_iter_max_LS=30; % Maximum number of iterations for Line Search
options.type_LS=2; % 1: Backtracking, 2: Matlab
options.TolX=1.e-2;     % For Matlab line search
options.alfa=0.3; % For backtracking
options.beta = .8; % For backtracking

%Setup the undeformed configuration
precompute;

history_E=[];
history_x=[];
history_delta=[];
history_F=[];

%loop on the load increments
for iload=1:length(lambda)
    %Define the boundary conditions
    x=x_eq;
    %if spring == 1
    %    force_sp=-K*abs((x(2*CC0'+1)-x(2*CC0'-1))).*x(2*CC0'-1)
    %    force(2*CC0'-1)=force_sp;
    %end
    load1.force = force*lambda(iload);
    switch example
        case {1, 2}
            x(1:2:end)=x_eq(1:2:end)*lambda(iload)/lambda(max(iload-1,1));
            load1.fixedvalues = x(load1.dofCC);
    end
    
    %x=x+rand(size(x))*.001; %random perturbations

    %Solve the equilibrium nonlinear system of equations
    [x_eq,iflag,iter,E_eq,Ensp] = Equilibrate(x,options,spring,K);
    [E_eq,grad_eq] = Energy(x_eq,3,Ensp);
    history_E(iload)=E_eq;
    history_x(iload,:)=x_eq;
    switch example
        case {0, 1, 2, 3}
            history_delta(iload)=x_eq(2*CC1(1)-1)-mesh1.x0(2*CC1(1)-1);
            history_F(iload)=sum(grad_eq(2*CC0'-1)); %Reaction
        case {4, 5}
            history_delta(iload)=mean(x_eq(dof_disp)-mesh1.x0(dof_disp));
            history_F(iload)=mean(load1.force(dof_force)); %Reaction
        otherwise
            disp('Case not implemented')
    end
    %Plot the deformed equilibrium configuration
    figure(1)
    clf
    DibujaMalla(mesh1.T,mesh1.x0,x_eq,'r',1)
    pause(0.01)
end


%Plot the deformation vs. force
figure(3)
plot(-(history_delta),abs(history_F),'ro-')
xlabel('\delta')
ylabel('Force')

%Solve the same problem using the linear theory of elasticity
[u,Reaction,delta] = linear_elasticity(codeLoad,lambda(end),example,CC1,dof_force,dof_disp,Ensp);
hold on
plot([0 -(delta)],[0,abs(Reaction)],'k-')
legend('Nonlinear elasticity','Linear elasticity')
figure(1)
hold on
DibujaMalla(mesh1.T,mesh1.x0,mesh1.x0+u,'k',1)
