function DibujaMalla(T,X_,XF_,str,nonum)
% DibujaMalla(T,X,str,nonum)
% X:    coordenadas nodales
% T:    conectividades (elementos)
% str:  s�mbolos y color para el dibujo (opcional)
% nonum = 1 si se quieren numeras los nodos (opcional)

X=[X_(1:2:end) X_(2:2:end)];
XF=[XF_(1:2:end) XF_(2:2:end)];

% definici�n del tipo de linea y color
if nargin == 2
  str1 = 'yo';
  str2 = 'y-';
else
  if str(1) == ':' | str(1) == '-'  
    str1 = 'yo'; 
    str2 = ['y' str];
  else
    str1 = [str(1) 'o'];
    str2 = str;
  end
end
 
nnodes = size(T,2);

% Dibujo de los elementos con la funci�n plot
if  size(X,2) == 2
  order = [1:nnodes,1];
  if nnodes == 8
    % definici�n del orden para un elemento de 8 nodos
    order = [1 5 2 6 3 7 4 8 1];
  end;
% Dibujo de los nodos
  plot(X(:,1),X(:,2),str1)
  hold on



% Dibujo de los elementos
  for j = 1:size(T,1)
    plot(X(T(j,order),1),X(T(j,order),2),'B')
    plot(XF(T(j,order),1),XF(T(j,order),2),str)
  end
end

% Numeraci�n de los nodos
  if nargin==4 
    if nonum==1
      for I=1:size(X,1)
        text(X(I,1),X(I,2),int2str(I))
      end
    end
  end

axis('equal')    
axis('off') 

hold off 
