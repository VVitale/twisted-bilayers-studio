function [supercell,psi_M1,psi_M2,psi_k1,psi_k2] = plot_exciton_wfc(index_orb_up,index_orb_down,r_h, orb_h_up, orb_h_down, Acvk, ...
         Cv, Cc, mcell, X, Y, Z, all_kpts, knum_tot, Nat, indc, indv, indk, ...
         dim, eig, outname, twisted)

% Define supercell in real space from k-grid
is = 0;
tot_r = knum_tot*Nat;
supercell = zeros(tot_r,3);
knum = sqrt(knum_tot);
psi_M1 = zeros(tot_r,1);
psi_M2 = zeros(tot_r,1);
nat_index = zeros(tot_r,1);
sum_M1 = 0.0;
sum_M2 = 0.0;
% tic
for in = 0 : knum - 1
   for jn = 0 : knum - 1
       R = in*mcell(1,:) + jn*mcell(2,:);
       for iat = 1 : Nat
           %is = (Nat+1)*in*(knum-1) + Nat*jn + iat;
           is = is + 1;
           r_e = [X(iat), Y(iat), Z(iat)];
           supercell(is,:) = r_e + R;
           nat_index(is) = iat;
       end
   end
end

parfor is = 1 : tot_r
      sum_M1 = 0.0;
      sum_M2 = 0.0;
      for idim = 1 : dim
             psi_vk_rh_up = conj(Cv([orb_h_up.Ham_index], indv(idim), indk(idim)));
             %psi_vk_rh_up = conj(Cv([index_orb_down{nat_index(is)}], indv(idim), indk(idim)));
             psi_ck_re_up = Cc([index_orb_up{nat_index(is)}], indc(idim), indk(idim));
             %psi_ck_re_up = Cc([orb_h_up.Ham_index], indc(idim), indk(idim));
             psi_vk_rh_down = conj(Cv([orb_h_down.Ham_index], indv(idim), indk(idim)));
             %psi_vk_rh_down = conj(Cv([index_orb_down{nat_index(is)}], indv(idim), indk(idim)));
             psi_ck_re_down = Cc([index_orb_down{nat_index(is)}], indc(idim), indk(idim));
             %psi_ck_re_down = Cc([orb_h_up.Ham_index], indc(idim), indk(idim));
             eikdotR = exp(1i*(dot(all_kpts(indk(idim),:), supercell(is,:) - r_h)));
             CvCc = (dot(psi_vk_rh_up,psi_vk_rh_up)*dot(psi_ck_re_up,psi_ck_re_up) + dot(psi_vk_rh_down,psi_vk_rh_down)*dot(psi_ck_re_down,psi_ck_re_down));
                       %sum(psi_vk_rh_up)*sum(psi_ck_re_down) + sum(psi_vk_rh_down)*sum(psi_ck_re_up));
             sum_M1 = sum_M1 + Acvk(idim,1) * eikdotR * CvCc;
             sum_M2 = sum_M2 + Acvk(idim,2) * eikdotR * CvCc;
      end
      psi_M1(is) = psi_M1(is) + abs(sum_M1)^2;
      psi_M2(is) = psi_M2(is) + abs(sum_M2)^2;
end
% toc
psi_M1 = psi_M1 / knum_tot^2;
psi_M2 = psi_M2 / knum_tot^2;

scatter3(supercell(:,1), supercell(:,2), supercell(:,3), 100, psi_M1, 'filled');
title(join(['$\epsilon =$', num2str(eig),' (eV)', '\textbf{r}$_h$ =', '(', num2str(r_h(1)),', ', num2str(r_h(2)),', ', num2str(r_h(3)),') \AA']),'Interpreter','latex','FontSize',22)
saveas(gcf, join([outname,'_M1.fig']));

scatter3(supercell(:,1), supercell(:,2), supercell(:,3), 100, psi_M2, 'filled');
title(join(['$\epsilon =$', num2str(eig),' (eV)', '\textbf{r}$_h$ =', '(', num2str(r_h(1)),', ', num2str(r_h(2)),', ', num2str(r_h(3)),') \AA']),'Interpreter','latex','FontSize',22)
saveas(gcf, join([outname,'_M2.fig']));

close all
if(~twisted)
   super_supercell = zeros(4*size(supercell,1),3);
   super_psi_M1 = zeros(4*size(psi_M1,1),1);
   super_psi_M2 = zeros(4*size(psi_M2,1),1);
   is = 0;
   for in = 0 : 1
      for jn = 0 : 1
         for inat = 1 : size(supercell,1)
            is = is + 1;
            super_supercell(is,:) = supercell(inat,:) + in*knum*mcell(1,:) + jn*knum*mcell(2,:);
            super_psi_M1(is,:) = psi_M1(inat);
            super_psi_M2(is,:) = psi_M2(inat);
         end
      end
   end
   scatter3(super_supercell(:,1),super_supercell(:,2),super_supercell(:,3),500,super_psi_M1,'filled')
   title(join(['$\epsilon =$', num2str(eig),' (eV)', '\textbf{r}$_h$ =', '(', num2str(r_h(1)),', ', num2str(r_h(2)),', ', num2str(r_h(3)),') \AA']),'Interpreter','latex','FontSize',22)
   saveas(gcf,join([outname,'_supercell_M1.fig']))
   scatter3(super_supercell(:,1),super_supercell(:,2),super_supercell(:,3),500,super_psi_M2,'filled')
   title(join(['$\epsilon =$', num2str(eig),' (eV)', '\textbf{r}$_h$ =', '(', num2str(r_h(1)),', ', num2str(r_h(2)),', ', num2str(r_h(3)),') \AA']),'Interpreter','latex','FontSize',22)
   saveas(gcf,join([outname,'_supercell_M2.fig']))
   close all
end
% Reciprocal space
psi_k1 = zeros(knum_tot,1);
psi_k2 = zeros(knum_tot,1);
for idim = 1 : dim
    psi_k1(indk(idim)) = psi_k1(indk(idim)) + abs(Acvk(idim,1)).^2;
    psi_k2(indk(idim)) = psi_k2(indk(idim)) + abs(Acvk(idim,2)).^2;
end
scatter3(all_kpts(:,1),all_kpts(:,2),zeros(size(psi_k1)),200,psi_k1(:),'filled')
saveas(gcf,join([outname,'_recspace_M1.fig']))
close all
scatter3(all_kpts(:,1),all_kpts(:,2),zeros(size(psi_k2)),200,psi_k2(:),'filled')
saveas(gcf,join([outname,'_recspace_M2.fig']))

clear psi_ck_rh psi_vk_rh  psi_ck_re psi_vk_re super_supercell super_psi_M1 super_psi_M2 CvCc eikdotR

end
