function [x_short t]= LineSearch(x_short,p,Ener,direction,options,spring,K)

t=1;
opts=optimset('TolX',options.TolX,'MaxIter',options.n_iter_max_LS);
t = fminbnd(@(t) Ener_1D(t,x_short,p,spring,K),0,2,opts);
x_short=x_short+t*p;