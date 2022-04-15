%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to plot Berry curvature and compute
% Chern number on a lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,F12] = chern(tb_vecs,knum_tot,phase1,phase2, ...
             multilayer,twisted,nspin,tot_norbs,noccs,all_kpts,b1,b2,outfname);
phi=pi/3;
if(multilayer)
   if(twisted)
      if(nspin==2)
         psi = complex(zeros(tot_norbs,2,knum_tot));
         for ik = 1 : knum_tot
            psi(:,1:2,ik) = [tb_vecs(:,noccs(ik),ik) tb_vecs(:,noccs(ik)-1,ik)];% tb_vecs(:,noccs(ik)-2,ik) tb_vecs(:,noccs(ik)-3,ik)]; % tb_vecs(:,noccs(ik)-4,ik) tb_vecs(:,noccs(ik)-5,ik) tb_vecs(:,noccs(ik)-6,ik) tb_vecs(:,noccs(ik)-7,ik)]; 
         end
      else
         psi = complex(zeros(tot_norbs,4,knum_tot));
         for ik = 1 : knum_tot
            psi(:,1:4,ik) = [tb_vecs(:,noccs(ik),ik) tb_vecs(:,noccs(ik)-1,ik) tb_vecs(:,noccs(ik)-2,ik) tb_vecs(:,noccs(ik)-3,ik)]; 
         end
      end
   else
      if(nspin==2)
         psi = complex(zeros(tot_norbs,4,knum_tot));
         for ik = 1 : knum_tot
            psi(:,1:4,ik) = [tb_vecs(:,noccs(ik),ik) tb_vecs(:,noccs(ik)-1,ik) tb_vecs(:,noccs(ik)-2,ik) tb_vecs(:,noccs(ik)-3,ik)]; 
         end
      else
         psi = complex(zeros(tot_norbs,2,knum_tot));
         for ik = 1 : knum_tot
            psi(:,1:2,ik) = [tb_vecs(:,noccs(ik),ik) tb_vecs(:,noccs(ik)-1,ik)]; 
         end
      end
   end
else
   if(nspin==2)
      psi = complex(zeros(tot_norbs,2,knum_tot));
      for ik = 1 : knum_tot
         psi(:,1:2,ik) = [tb_vecs(:,noccs(ik),ik) tb_vecs(:,noccs(ik)-1,ik)]; 
      end
   else
      psi = complex(zeros(tot_norbs,1,knum_tot));
      for ik = 1 : knum_tot
         psi(:,1,ik) = tb_vecs(:,noccs(ik),ik);
      end
   end
end
F12 = nabelian_berry_curv(psi,knum_tot,phase1,phase2);
% Check  -pi > F12 >= pi 
ind = find(abs(F12) >= pi);
if(any(ind))
   ind
   error('The absolute value of Berry curvature is greater than pi. Aborting ...')
end
% Filter out small values
ind = find(abs(F12) < pi*0.00001);
for i = 1 : knum_tot
   if(any(ismember(ind,i)))
      F12(i) = complex(0.0);
   else
      F12(i) = F12(i);
   end
end
close all
dk = norm(b1)/sqrt(knum_tot);
F12 = F12/dk^2;
F12 = -1i*F12;
supercell = [all_kpts; all_kpts + b1; all_kpts + b2; all_kpts + b1 + b2;
             all_kpts + 2*b1; all_kpts + 2*b2; all_kpts + 2*b1 + b2;
             all_kpts + b1 + 2*b2; all_kpts + 2*(b1 + b2);
             all_kpts + 3*b1; all_kpts + 3*b2; all_kpts + 3*b1 + b2;
             all_kpts + 3*b2 + b1; all_kpts + 3*b1 + 2*b2;
             all_kpts + 3*b2 + 2*b1; all_kpts + 3*(b1+b2)];
sF12 = [real(F12); real(F12); real(F12); real(F12);
        real(F12); real(F12); real(F12); real(F12);
        real(F12); real(F12); real(F12); real(F12);
        real(F12); real(F12); real(F12); real(F12)];
scatter(supercell(:,1),supercell(:,2),100,sF12,'filled')
saveas(gca,join([outfname,'_Berry_VB.fig']))
close all
centre = 2*(b1+b2);
x = linspace(centre(1)-norm((b1+b2)/3)-0.05,centre(1)+norm((b1+b2)/3)+0.05,300); 
y = linspace(centre(2)-0.5*(b2(2))-0.05,centre(2)+0.5*(b2(2))+0.05,300);
[X,Y] = meshgrid(x,y);
f = scatteredInterpolant(supercell(:,1),supercell(:,2),sF12,'natural');
Z = f(X,Y);
surf(X,Y,Z)
shading interp
%cmap = cbrewer('seq','YlGnBu',20,'PCHIP');
colormap(jet);
colorbar;
axis([x(1) x(end) y(1) y(end)])
hold on
% Draw exagon
Rotphi = [cos(pi/3) -sin(pi/3) 0;sin(pi/3) cos(pi/3) 0; 0 0 1];
k1=centre+(2*b1-b2)/3;
k2 = (k1-centre)*Rotphi' + centre;
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
%set(gca,'XtickLabel',[-(2*b1(1)-b1(2))/3,0.0,(2*b1(1)-b1(2))/3],'FontSize',15);
%set(gca,'YtickLabel',[-b1(1)/2,0.0,b1(1)/2],'FontSize',15);
height=1600;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
text(centre(1),centre(2),height, '\Gamma','FontSize',25,'Color','w')
text(k1(1),k1(2),height, 'K_{+}','FontSize',25,'Color','w')
text(k2(1),k2(2),height, 'K_{-}','FontSize',25,'Color','w')
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
k1 = k2;
k2 = (k1-centre)*Rotphi' + centre;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
k1 = k2;
k2 = (k1-centre)*Rotphi' + centre;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
k1 = k2;
k2 = (k1-centre)*Rotphi' + centre;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
k1 = k2;
k2 = (k1-centre)*Rotphi' + centre;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
k1 = k2;
k2 = (k1-centre)*Rotphi' + centre;
x = linspace(k1(1),k2(1),height);
y = linspace(k1(2),k2(2),height);
plot3(x,y,height*ones(size(x)),'w','LineWidth',3)
xlabel('k_x','FontSize',50);
ylabel('k_y','FontSize',50);
saveas(gca,join([outfname,'_Berry_VB_interp.fig']))
% Compute Chern number
c = 1/(2*pi)*sum(F12,1)*dk^2;
% Check c is an integer
if(mod(fix(real(c)),1) ~= 0)
   error('Chern number is not an integer!')
else
   c = fix(real(c))
end
clear ind phase1 phase2 psi;

% Compute Berry curvature
function F12 = nabelian_berry_curv(unk,tot_knum,phase1,phase2)

Nk1 = sqrt(tot_knum);
F12 = complex(zeros(tot_knum,1));
% Run over direction 2 first
% All k-points but last Nk1
for kl = 1 : tot_knum - Nk1
   if(mod(kl,Nk1) == 0)
      U2  = det(unk(:,:,kl)'*(phase2.*unk(:,:,kl-Nk1+1)));
      U1  = det(unk(:,:,kl)'*unk(:,:,kl+Nk1));
      U12 = det((phase2.*unk(:,:,kl-Nk1+1))'*(phase2.*unk(:,:,kl+1)));
      U21 = det(unk(:,:,kl+Nk1)'*(phase2.*unk(:,:,kl+1)));
   else
      U2  = det(unk(:,:,kl)'*unk(:,:,kl+1));
      U1  = det(unk(:,:,kl)'*unk(:,:,kl+Nk1));
      U12 = det(unk(:,:,kl+1)'*unk(:,:,kl+Nk1+1));
      U21 = det(unk(:,:,kl+Nk1)'*unk(:,:,kl+1+Nk1));
   end
   U1 = U1/norm(U1);
   U2 = U2/norm(U2);
   U21 = U21/norm(U21);
   U12 = U12/norm(U12);
   F12(kl) = log(U1*U21*inv(U12)*inv(U2));
end
jkl = 0;
for kl = tot_knum - Nk1 + 1 : tot_knum
   jkl = jkl + 1;
   if(mod(kl,Nk1)==0)
      U2  = det(unk(:,:,kl)'*(phase2.*unk(:,:,kl-Nk1+1)));
      U1  = det(unk(:,:,kl)'*(phase1.*unk(:,:,jkl)));
      U12 = det((phase2.*unk(:,:,kl-Nk1+1))'*(phase2.*phase1.*unk(:,:,jkl-Nk1+1)));
      U21 = det((phase1.*unk(:,:,jkl))'*(phase2.*phase1.*unk(:,:,jkl-Nk1+1)));
   else
      U2  = det(unk(:,:,kl)'*unk(:,:,kl+1));
      U1  = det(unk(:,:,kl)'*(phase1.*unk(:,:,jkl)));
      U12 = det(unk(:,:,kl+1)'*(phase1.*unk(:,:,jkl+1)));
      U21 = det((phase1.*unk(:,:,jkl))'*(phase1.*unk(:,:,jkl+1)));
   end
   U1 = U1/norm(U1);
   U2 = U2/norm(U2);
   U21 = U21/norm(U21);
   U12 = U12/norm(U12);
   F12(kl) = log(U1*U21*inv(U12)*inv(U2));
end
end
% function F12 = abelian_berry_curv(unk,tot_knum,phase1,phase2)
% 
% Nk1 = sqrt(tot_knum);
% F12 = complex(zeros(tot_knum,1));
% % Run over direction 2 first
% % All k-points but last Nk1
% for kl = 1 : tot_knum - Nk1
%    if(mod(kl,Nk1) == 0)
%       U2  = dot(unk(:,kl),(phase2.*unk(:,kl-Nk1+1)));
%       U1  = dot(unk(:,kl),unk(:,kl+Nk1));
%       U12 = dot((phase2.*unk(:,kl-Nk1+1)),(phase2.*unk(:,kl+1)));
%       U21 = dot(unk(:,kl+Nk1),(phase2.*unk(:,kl+1)));
%    else
%       U2  = dot(unk(:,kl),unk(:,kl+1));
%       U1  = dot(unk(:,kl),unk(:,kl+Nk1));
%       U12 = dot(unk(:,kl+1),unk(:,kl+Nk1+1));
%       U21 = dot(unk(:,kl+Nk1),unk(:,kl+1+Nk1));
%    end
%    U1 = U1/norm(U1);
%    U2 = U2/norm(U2);
%    U21 = U21/norm(U21);
%    U12 = U12/norm(U12);
%    F12(kl) = log(U1*U21*inv(U12)*inv(U2));
% end
% jkl = 0;
% for kl = tot_knum - Nk1 + 1 : tot_knum
%    jkl = jkl + 1;
%    if(mod(kl,Nk1)==0)
%       U2  = dot(unk(:,kl),(phase2.*unk(:,kl-Nk1+1)));
%       U1  = dot(unk(:,kl),(phase1.*unk(:,jkl)));
%       U12 = dot((phase2.*unk(:,kl-Nk1+1)),(phase2.*phase1.*unk(:,jkl-Nk1+1)));
%       U21 = dot((phase1.*unk(:,jkl)),(phase2.*phase1.*unk(:,jkl-Nk1+1)));
%    else
%       U2  = dot(unk(:,kl),unk(:,kl+1));
%       U1  = dot(unk(:,kl),(phase1.*unk(:,jkl)));
%       U12 = dot(unk(:,kl+1),(phase1.*unk(:,jkl+1)));
%       U21 = dot((phase1.*unk(:,jkl)),(phase1.*unk(:,jkl+1)));
%    end
%    U1 = U1/norm(U1);
%    U2 = U2/norm(U2);
%    U21 = U21/norm(U21);
%    U12 = U12/norm(U12);
%    F12(kl) = log(U1*U21*inv(U12)*inv(U2));
% end
% end
end
