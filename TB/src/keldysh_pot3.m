function y = keldysh_pot2(q,R,positions,nat,a,shift,knum_tot)
% These are hard-coded at the moment
% TODO: make these input variables
ed = 2.5;
r0 = 33.875/ed;
U = 1;
keldysh_prefactor = 22.59;
y = zeros(nat);
dist_mat = zeros(nat);
Tx = [positions(:).x] - shift(1);
Ty = [positions(:).y] - shift(2);
Tz = [positions(:).z] - shift(3);
Tx_mat = repmat(Tx',[nat,1]) - repmat(Tx,[1,nat]);
Ty_mat = repmat(Ty',[nat,1]) - repmat(Ty,[1,nat]);
Tz_mat = repmat(Tz',[nat,1]) - repmat(Tz,[1,nat]);
for iR = 1 : knum_tot
   eiqdotR = exp(1i*dot(q,R(iR,:)));
   RTx_mat = (Tx_mat + R(iR,1)).^2;
   RTy_mat = (Ty_mat + R(iR,2)).^2;
   RTz_mat = (Tz_mat + R(iR,3)).^2;
   dist_mat = sqrt(RTx_mat + RTy_mat + RTz_mat);
   dist_mat(dist_mat == 0) = a/sqrt(3);
   eiqdotT = exp(1i.*(q(1).*Tx_mat + q(2).*Ty_mat + q(3).*Tz_mat));
   tV  = -U/(ed*r0) * (StruveH0(dist_mat/r0)-bessely(0,dist_mat/r0));
   y = keldysh_prefactor*eiqdotR.*eiqdotT.*tV;
end
end
