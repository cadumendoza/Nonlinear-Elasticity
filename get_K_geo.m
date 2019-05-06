function [Hess_E] = get_K_geo(x)
global mod1 mesh1 load1 el1 undeformed1

% We: deformation energy
% X: initial nodal coordinates matrix
% XF: final nodal coordinates matrix
% T: connectivities matrix
% pospg: Gauss points coordinates (reference elements)
% pespg: Gauss points weights (reference elements)
% N: shape functions (reference elements) on Gauss points
% Nxi, Neta: derivatives of N with respect to xi and eta on Gauss points

nelem=size(mesh1.T,1);
nnod =size(mesh1.T,2);
npoin = size(x,1)/2;
Hess_E = [];
ig=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
Hess_E = zeros(npoin*2,npoin*2);
for ielem = 1:nelem
    a=2*mesh1.T(ielem,:)-1;
    b=2*mesh1.T(ielem,:);
    Te = [ a(1) b(1) a(2) b(2) a(3) b(3) a(4) b(4) ];
    xe = x(Te);
    grad_Ee = zeros(1,8);
    Hess_Ee = zeros(8);
    for igaus = 1:el1.ngaus
        ig=ig+1;
        defGrad=[xe(1:2:end)' ; xe(2:2:end)']*undeformed1.DNDX{ig}';
        C=defGrad'*defGrad;
        [W,S,CC]=NeoHookean_3([C(1,1) C(2,2) C(1,2)],mod1.lambda,mod1.mu,3);
        DNDX=undeformed1.DNDX{ig};            
        temp=[defGrad(:,1)*DNDX(1,:) defGrad(:,2)*DNDX(2,:) ...
              defGrad(:,2)*DNDX(1,:)+defGrad(:,1)*DNDX(2,:)];
        B0=reshape(temp,8,3)';
        S_=[S(1) S(3); S(3) S(2)];
        H=DNDX'*S_*DNDX;
        HH=zeros(8);
        HH(1:2:8,1:2:8)=H;
        HH(2:2:8,2:2:8)=H;
        Hess_Ee = Hess_Ee + (HH) * undeformed1.dvol0(ig);              
    end
   Hess_E(Te,Te) = Hess_E(Te,Te) + Hess_Ee;
end



%%%%%%%%%%%%%%%%%%%%


