function [X,T]=CreaMalla(x1,x2,y1,y2,nx,ny)
% [X,T]=CreaMalla(x1,x2,y1,y2,nx,ny)
% Crea la topologia de una malla estructurada de cuadril�teros
% en un dominio rectangular [x1,x2]x[y1,y2]
% nx,ny: n�mero de divisiones en cada direcci�n


X = zeros((nx+1)*(ny+1),2);
T = zeros(nx*ny,4);

hx = (x2-x1)/nx;
hy = (y2-y1)/ny;

xs = [x1:hx:x2]';
unos = ones(nx+1,1);

% Coordenadas de los nodos
for i=1:ny+1
 ys = ((i-1)*hy+y1)*unos;
 posi = [(i-1)*(nx+1)+1:i*(nx+1)]; 
 X(posi,:)=[xs,ys];
end


% Conectividades
for i=1:ny
 for j=1:nx
   ielem = (i-1)*nx+j;
   inode = (i-1)*(nx+1)+j;
   T(ielem,:) = [ inode inode+1 inode+(nx+2) inode+(nx+1)];
 end
end
