function y = keldysh_pot2(q,R,positions,nat,a,shift,knum_tot)
% These are hard-coded at the moment
% TODO: make these input variables
ed = 2.5;
r0 = 33.875/ed;
U = 1;
keldysh_prefactor = 22.59;
y = zeros(nat);
qmat = repmat(q,[nat,1]);
T = [positions(:).x,positions(:).y,positions(:).z] - shift;
tij_mat = repmat(T,[1,nat]) - repmat(reshape(T',[1,nat*3]),[nat,1]);
for iR = 1 : knum_tot
   %qdotR = dot(q,R(iR,:));
   for inat = 1 : nat
      %ti = [positions.x(inat),positions.y(inat),positions.z(inat)] - shift;
      tijR_mat = tij_mat(:,(inat-1)*3 + 1 : (inat-1)*3 + 3 )  + R(iR,:);
      eiqdotRt_mat = exp(1i.*dot(qmat,tijR_mat,2));
      dist_mat = vecnorm(tijR_mat,2,2);
      dist_mat(dist_mat == 0) = a/sqrt(3);
     % for jnat = 1 : nat
     %     tj = [positions.x(jnat),positions.y(jnat),positions.z(jnat)] - shift;
     %     tji = tj - ti;
     %     qdott = dot(q,tji);
     %     dist = norm(R(iR,:) + tji);
     %     tV = 0;
     %     if(dist==0)
     %        tV  = -U/(ed*r0) * (StruveH0(a/sqrt(3)/r0)-bessely(0,a/sqrt(3)/r0));
     %     else
             tV  = -U/(ed*r0) * (StruveH0(dist_mat/r0)-bessely(0,dist_mat/r0));
     %     end
     %     y(inat,jnat) = keldysh_prefactor*tV*exp(1i*qdotR)*exp(1i*qdott);
          y(inat,:) = keldysh_prefactor.*tV'.*conj(eiqdotRt_mat');
     % end
   end
end
end
