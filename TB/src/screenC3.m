function u = screenC3(q,Nk,nat,structure,mcell,a)
u = zeros(nat);
M = complex(zeros(nat));
shift = 2/3*mcell(1,:)+1/3*mcell(2,:);
Ni = -(Nk-1)/2;
Nf = -Ni;
dN = 1;
kn = 0;
for im = Ni : dN : Nf
    for jn = Ni : dN : Nf
        kn = kn + 1;
        R = im*mcell(1,:) + jn*mcell(2,:);
        VR_ij = zeros(nat);
        M = zeros(nat);
        for inat = 1 : nat
            ti = [structure.x(inat),structure.y(inat),structure.z(inat)] - shift;
            for jnat = 1 : nat
                tj = [structure.x(jnat),structure.y(jnat),structure.z(jnat)] - shift;
                qdotRij = dot(q,tj - ti);
                M(inat,jnat) = exp(1i*qdotRij);
                dist = norm(R + tj - ti);
                VR_ij(inat,jnat) = keldysh_pot(dist,inat,jnat,a);
            end
        end
        u = u + exp(1i*dot(q,R)).*M.*VR_ij;
   end
end
clear R qdotR
end


