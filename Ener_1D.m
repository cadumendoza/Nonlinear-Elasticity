function energy = Ener_1D(t,x,p,spring,K,iload)

energy = Ener_short(x+t*p,1,spring,K,iload);
