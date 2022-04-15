function y = keldysh_pot4(Rnorm,a)
ed = 2.5;
r0 = 33.875/ed;
U = 1;
keldysh_prefactor = 22.59;
tV = 0;
if(Rnorm==0)
   tV  = -U/(ed*r0) * (StruveH0(a/sqrt(3)/r0)-bessely(0,a/sqrt(3)/r0));
else
   tV  = -1/(ed*r0) * (StruveH0(Rnorm/r0)-bessely(0,Rnorm/r0));
end
y = keldysh_prefactor*tV;
end
