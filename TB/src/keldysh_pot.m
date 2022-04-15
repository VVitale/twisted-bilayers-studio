function y = keldysh_pot(dist,inat,jnat,a)
ed = 2.5;
r0 = 33.875/ed;
U = 1;
keldysh_prefactor = 22.59;
tV = 0;
if(dist==0)
   tV  = -U/(ed*r0) * (StruveH0(a/sqrt(3)/r0)-bessely(0,a/sqrt(3)/r0));
else
   tV  = -1/(ed*r0) * (StruveH0(dist/r0)-bessely(0,dist/r0));
end
y = keldysh_prefactor*tV;
end
