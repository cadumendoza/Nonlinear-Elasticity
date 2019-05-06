function [pospg,pespg]=Cuadratura (ngaus) 
% [pospg,pespg]=Cuadratura (ngaus) 
% pospg, pespg: posici�n y pesos de los puntos de Gauss en el 
%               elemento de referencia 
% ngaus:        n�mero de puntos de Gauss del elemento 
 
if (ngaus==4) 
   pos1 = 1/sqrt(3); 
   pospg=[-pos1 -pos1 
           pos1 -pos1 
           pos1  pos1 
          -pos1  pos1]; 
   pespg=[ 1 1 1 1]; 
else 
   error(' Cuadratura no contemplada') 
end 
 
           
