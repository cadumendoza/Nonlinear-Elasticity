function precompute
global mesh1 el1 undeformed1

nelem=size(mesh1.T,1);

ig=0;
for ielem = 1:nelem
    a=2*mesh1.T(ielem,:)-1;
    b=2*mesh1.T(ielem,:);
    Te = [ a(1) b(1) a(2) b(2) a(3) b(3) a(4) b(4) ];
    Xe = mesh1.x0(Te);
    for igaus = 1:el1.ngaus
        ig=ig+1;
        N_igaus = el1.N(igaus,:);
        Nxi_igaus = el1.DN(igaus,:,1);
        Neta_igaus = el1.DN(igaus,:,2);
        DN(1,:)=el1.DN(igaus,:,1);
        DN(2,:)=el1.DN(igaus,:,2);
        dXdxi  = [Nxi_igaus*(Xe(1:2:end))	Neta_igaus*(Xe(1:2:end))
                  Nxi_igaus*(Xe(2:2:end))	Neta_igaus*(Xe(2:2:end))];
        temp=inv(dXdxi);
        undeformed1.dxidX{ig}=temp;
        undeformed1.dvol0(ig)=el1.pespg(igaus)*det(dXdxi);
        undeformed1.DNDX{ig}=undeformed1.dxidX{ig}'*DN;
    end
end
