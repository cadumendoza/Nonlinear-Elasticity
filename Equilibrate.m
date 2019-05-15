function [x_equil,iflag,iter,Ener] = Equilibrate(x,options,spring,K)
global mod1 mesh1 load1 el1 undeformed1


err_plot=[];
err_plot1=[];
[x_short] = short(x);

switch options.method
  case 0,   %vanilla Newton-Raphson
    iter=0;
    err_x=100;
    err_f=100;
    [Ener,grad_E,Hess_E] = Ener_short(x_short,3,spring,K);
    while (iter<=options.n_iter_max) & ...
        ( (err_x>options.tol_x) | ...
        (err_f>options.tol_f))
      iter=iter+1;
      dx = -Hess_E\grad_E;
      x_short=x_short+dx;
      [Ener,grad_E,Hess_E] = Ener_short(x_short,3,spring,K);
      err_x=norm(dx)/norm(x_short);
      err_f=norm(grad_E);
      err_plot=[err_plot err_x];
      err_plot1=[err_plot1 err_f];
      %fprintf('Iteration %i, errors %e %e \n', iter,err_x,err_f)
    end
    %Check positive definiteness
    if options.info==3
      [V,D] = eig(Hess_E);
      D=diag(D);
      if ((min(D))<=-1e-6*abs(max(D)))
           fprintf('Warning, the Hessian has a negative eigenvalue \n')
      end
    end
  case 1,   % Newton-Raphson with line search
    iter=0;
    err_x=100;
    err_f=100;
    [Ener,grad_E,Hess_E] = Ener_short(x_short,3,spring,K);
    while (iter<=options.n_iter_max) & ...
        ( (err_x>options.tol_x) | ...
        (err_f>options.tol_f))
      iter=iter+1;
      ihess=1;
      for Ihess=2:2:length(load1.dofSp)
          Hess_E(Ihess,Ihess)=Hess_E(Ihess,Ihess)+load1.Ks(ihess);
          ihess=ihess+1;
      end
      for Ihess=2:2:length(load1.dofSp)*2
          grad_E(Ihess)=grad_E(Ihess)-load1.forcesp(Ihess);
      end
      dx = -Hess_E\grad_E;        
    
      if options.linesearch==1
          direction=dx'*grad_E;
          if (direction)>0 
              disp('reverse direction')
              dx=-dx;
              direction=-direction;
          end
          [x_short, t]=LineSearch(x_short,dx,Ener,direction,options,spring,K);
      else
          t=1;
          x_short=x_short+dx;
      end
      
      [Ener,grad_E,Hess_E] = Ener_short(x_short,3,spring,K);
      err_x=abs(t)*norm(dx)/norm(x_short);
      err_f=norm(grad_E);
      err_plot=[err_plot err_x];
      err_plot1=[err_plot1 err_f];
      %fprintf('Iteration %i, errors %e %e \n', iter,err_x,err_f)
    end
    %Check positive definiteness
    if options.info==3
      [V,D] = eig(Hess_E);
      D=diag(D);
      if ((min(D))<=-1e-6*abs(max(D)))
           fprintf('Warning, the Hessian has a negative eigenvalue \n')
      end
    end
    otherwise,
    error('This option does not exist');
end

[x_equil] = long(x_short);
if (iter<=options.n_iter_max)
  if options.info>0
    fprintf('The equilibration was successful in %i iterations \n',iter)
  end
  iflag=1;
else
  if options.info>0
    fprintf('The equilibration was not reached in %i iterations \n',iter)
  end
  iflag=0;
end

%Output
if options.info>=2
  figure(5)
  hold off
  semilogy([1:iter],err_plot,'ro-',[1:iter],err_plot1,'bo-'); % rojo desplazamiento, azul fuerza
  xlabel('no. iterations')
  ylabel('log(Error)')
  legend('Displacement','Force')
  %pause
end



