function u = screenC2(q,Nk,nat,VR,lattice,structure,mcell)
u = zeros(nat);
M = complex(zeros(nat));
shift = 2/3*mcell(1,:)+1/3*mcell(2,:);
for kn = 1 : Nk
    R = lattice(kn,:);
    for inat = 1 : nat
        ti = [structure.x(inat),structure.y(inat),structure.z(inat)] - shift;
        for jnat = 1 : nat
            tj = [structure.x(jnat),structure.y(jnat),structure.z(jnat)] - shift;
            qdotRij = dot(q,tj - ti);
            M(inat,jnat) = exp(1i*qdotRij);
        end
    end
    u = u + exp(1i*dot(q,R)).*M.*VR{kn,1};
end
clear R qdotR
end
