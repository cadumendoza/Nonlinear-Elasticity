function [N,Nxi,Neta]=FuncForm(pospg) 
% [N,Nxi,Neta]=FuncForm(pospg) 
% N, Nxi, Neta: matrices que almacenan los valores de las funciones de forma 
%               y derivadas (en coordenadas locales del elemento de referencia)
%               en los puntos de Gauss (cada punto de Gauss en una fila) 
% pospg:        posición de los puntos de Gauss en el elemento de referencia 
 
xi = pospg(:,1); eta = pospg(:,2); 
N    = [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4]; 
Nxi  = [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4]; 
Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ]; 
