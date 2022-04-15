function u = screenC4(q,Nk,mcell,a)
Ni = -(Nk-1)/2;
Nf = -Ni;
dN = 1;
kn = 0;
u = 0.0;
for im = Ni : dN : Nf
    for jn = Ni : dN : Nf
        kn = kn + 1;
        R = im*mcell(1,:) + jn*mcell(2,:);
        eiqdotR = exp(1i*dot(q,R));
        Rnorm = norm(R);
        u = u + eiqdotR*keldysh_pot4(Rnorm,a);
   end
end
clear R qdotR
end


