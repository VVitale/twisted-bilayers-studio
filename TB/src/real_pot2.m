function [VR,lattice] = real_pot2(mcell,structure,a,nat,Nk)
% Epsilon, r0 and U (these are hard-coded at the moment)
ed = 2.5;
r0 = 33.875/ed;
U = 1;
keldysh_prefactor = 22.59;
lattice=zeros(Nk*Nk,3);
VR = cell(Nk,1);
kn = 0;
if(mod(Nk,2)==0)
    Ni = -(Nk-1)/2;
    Nf = -Ni;
    dN = 1;
else
    Ni = -(Nk-1)/2;
    Nf = -Ni;
    dN = 1;
end

% Shift centre of unit cell
shift = 2/3*mcell(1,:)+1/3*mcell(2,:);

% Compute the keldysh potential VR in real space
for im = Ni : dN : Nf
    for jn = Ni : dN : Nf
        kn = kn + 1;
        lattice(kn,:) = im*mcell(1,:) + jn*mcell(2,:);
        tV = zeros(nat,nat);
        for inat = 1 : nat
            ti = [structure.x(inat),structure.y(inat),structure.z(inat)] - shift;
            for jnat = 1 : nat
                tj = [structure.x(jnat),structure.y(jnat),structure.z(jnat)] - shift;
                dist = norm(lattice(kn,:) + tj - ti);
                if(dist==0)
                   tV(inat,jnat)  = -U/(ed*r0) * (StruveH0(a/sqrt(3)/r0)-bessely(0,a/sqrt(3)/r0));
                else
                   tV(inat,jnat)  = -1/(ed*r0) * (StruveH0(dist/r0)-bessely(0,dist/r0));
                end
                
            end
        end
        VR{kn} = keldysh_prefactor.*tV;
    end
end
clear R dist
end
