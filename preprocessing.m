function [dof_force, dof_disp, lambda, x_eq, CC0, CC1, force, codeLoad]=preprocessing(example,material)
global mod1 mesh1 load1 el1 undeformed1



dof_force=[];
dof_disp=[];
switch example
    case 0 %upsetting of a block, dead load
        x1=0;x2=1;
        y1=0;y2=1;
        nx=10;ny=10;
        lambda=[0:0.03:1.]; codeLoad=0;  % 0 for applied force and 1 for applied displacement
        codeLoad=0;
        mod1.force = 5e0;
    case 1 %upsetting of a block, imposed displacements
        x1=0;x2=1;
        y1=0;y2=1;
        nx=10;ny=10;
        lambda=[1:.025:2];
        %lambda=[1:-.01:0.5];
        codeLoad=1;
        mod1.force = 0;
    case 2 %compression of a slender beam, imposed displacements
        x1=0;x2=40;
        y1=0;y2=1;
        nx=40;ny=3;
        lambda=[1:-.005:0.5];
        codeLoad=1;
        mod1.force = 0;
    case 3 %compression of a slender beam, dead load
        x1=0;x2=40;
        y1=0;y2=1;
        nx=40;ny=3;
        lambda=-[0:0.002:.024];
        codeLoad=0;
        mod1.force = 100;
    case 4 % arch, dead load at center of the arch
        x1=0;x2=40;
        y1=0;y2=1;
        nx=40;ny=3;
        codeLoad=1;
        Rad=30;
        dof_force=2*(nx/2+1); lambda=-[0:0.004:.048];
        dof_disp=2*(nx/2+1);
        mod1.force = 1;
    case 5 % arch, dead load near the supports
        x1=0;x2=40;
        y1=0;y2=1;
        nx=40;ny=3;
        codeLoad=1;
        Rad=30;
        dof_force=[2*(nx/2+1-16) 2*(nx/2+1+16) ]; lambda=-[0:0.004:.048]*12;
        dof_disp=2*(nx/2+1);
        mod1.force = 1;
    otherwise
        disp('Case not implemented')
end

%Make the mesh
[XA,mesh1.T]=CreaMalla(x1,x2,y1,y2,nx,ny);
mesh1.x0=zeros(2*size(XA,1),1);
mesh1.x0(1:2:end) =  XA(:,1);
mesh1.x0(2:2:end) =  XA(:,2);


switch example
    case {4, 5}
        x_=mesh1.x0;
        theta=(x_(1:2:end)-20)/(Rad);
        x_(2:2:end) = (Rad+mesh1.x0(2:2:end)).*(cos(theta));
        x_(1:2:end) = (Rad+mesh1.x0(2:2:end)).*sin(theta);
        mesh1.x0=x_;
end


%Element
el1.ngaus=4;
[el1.pospg,el1.pespg]=Cuadratura (el1.ngaus);
[el1.N,el1.Nxi,el1.Neta]=FuncForm(el1.pospg);
el1.DN(:,:,1)=el1.Nxi;
el1.DN(:,:,2)=el1.Neta;

%Model parameters
switch material
    case 1
      mod1.potential = 1; % 1 is Neohookean
      mod1.mu=1;
      mod1.lambda=100; % not too large
    case 2
      mod1.potential = 2; % 2 is transversely isotropic model
      mod1.mu=1.;
      mod1.c0=80;
      mod1.c1=5;
      mod1.kappa=100;
      theta=pi/4;%pi/4;
      mod1.N_fib=[cos(theta); sin(theta)];
    case 3
      mod1.potential = 3; % 3 is KsV
      mod1.mu=1;
      mod1.lambda=100; % not too large
    otherwise
        error('Material not implemented')
end

%Initial Configuration
x_eq=mesh1.x0;

%Loading
%Fixed dofs
CC0=1:(nx+1):((ny)*(nx+1)+1);
CC1=(nx+1):(nx+1):((ny+1)*(nx+1));
dofCC0=[2*CC0'-1
    2*CC0'];
if codeLoad==0          % dead load, confined
    load1.dofCC=[dofCC0; 2*CC1'];
    dofCC1=[2*CC1'-1];
elseif codeLoad==1      % imposed displacement, confined
    dofCC1=[2*CC1'-1
        2*CC1'];
    load1.dofCC=[dofCC0;dofCC1];
elseif codeLoad==2      % dead load, linear springs
    
elseif codeLoad==3      % dead load, linear springs
    
else
    error('Loading not impemented')
end
load1.fixedvalues = x_eq(load1.dofCC);
load1.ndofCC=size(load1.dofCC,1);
load1.dofFree=[1:length(mesh1.x0)];
load1.dofFree(load1.dofCC) = [];
% Forces
force = zeros(size(mesh1.x0));
switch example
    case {0, 3}
        force_nod=mod1.force/ny;
        force(dofCC1)=force_nod;
        force(dofCC1(1))=force(dofCC1(1))/2;
        force(dofCC1(end))=force(dofCC1(end))/2;
    case {4, 5}
        force(dof_force)=mod1.force;
end

