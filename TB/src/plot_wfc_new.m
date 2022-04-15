%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to plot wave function on atoms for
% a given k-point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% tb_vecs - TB eigenvectors
% orbitals - orbitals class
% X, Y, Z  - Cartesian coordinates of atoms
% at_name - ordered list of atomic symbols
% mcell - lattice vectors of moire cell (ang)
% natoms - Number of atoms
% norbs - Number of orbitals
% n - band index of psi state
% kpoint - k-point in MP grid
% filename - Output file name
%
% Outputs:
% psink - wave function at k squared
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psink= plot_wfc(tb_vecs,orbitals,X,Y,Z,at_name,mcell,natoms,norbs,filename,twisted,multilayer,nlayer,k)

psink = zeros(natoms,1);
for ir = 1 : natoms
   r = [X(ir),Y(ir),Z(ir)];
   orb1 = findobj(orbitals,'Centre',r);
   orb2 = findobj(orbitals,'Pcentre',r);
   if(~isempty(orb1))
      iorb = orb1; %orb1.Ham_index;
      %centre = orbitals(iorb).Centre;
   elseif(~isempty(orb2))
      iorb = orb2; %orb2.Ham_index;
      %centre = orbitals(iorb).Pcentre;
   end
   if(twisted && multilayer)
       %psink(ir) = psink(ir) + dot(tb_vecs([iorb.Ham_index]),tb_vecs([iorb.Ham_index]));
       psink(ir) = psink(ir) + abs(sum(tb_vecs([iorb.Ham_index])))^2;
   end
end
clear orb1 orb2 iorb
nsublayers = 3*nlayer;
llim = 0;
ulim = 2.5;
for isl = 1 : nsublayers
   ilayer(isl,:) = find(Z > llim & Z < ulim);
   llim = llim + 2;
   ulim = ulim + 2;
end
pos = cell(nsublayers,1);
for il = 1 : nsublayers
   pos{il} = [X(ilayer(il,:)),Y(ilayer(il,:)),Z(ilayer(il,:))];
end
mean_z = mean(Z);
figure(1)
colormap(hot);
supercell = cell(nsublayers,1); %zeros(6,9*size(psink,1),3);
spsink = cell(nsublayers,1); %zeros(6,9*size(psink,1),1);
for il = 1 : nsublayers
   jp = 0;
   tmp = zeros(size(ilayer(il,:),2),3);
   tmp2 = zeros(size(ilayer(il,:),2),1);
   tmp3 = zeros(size(ilayer(il,:),2),3);
   for n = - 2 : 2
      for m = -2 : 2
          for ip = 1 : size(ilayer(il,:),2)
             jp = jp + 1;
             tmp3 = pos{il};
             tmp(jp,:) = tmp3(ip,:) + n*mcell(1,:) + m*mcell(2,:);
             tmp2(jp) = psink(ilayer(il,ip))/(tmp(jp,3)-mean_z)^2;
          end
      end
end
   supercell{il} = tmp;
   spsink{il} = tmp2;
end
clear tmp tmp2 tmp3;
phi = acos(mcell(1,1)/norm(mcell(1,:)));
Rphi = [cos(-phi) -sin(-phi) 0;sin(-phi) cos(-phi) 0; 0 0 1];
rot_cell = mcell*Rphi';
x = linspace(-rot_cell(1,1),rot_cell(1,1),300); 
y = linspace(-rot_cell(2,2),rot_cell(2,2),300); 
[XX,YY] = meshgrid(x,y);
ZZ = 0.*XX + 0.*YY;
for il = 1 : nsublayers
spos = supercell{il}*Rphi';
f = scatteredInterpolant(spos(:,1),spos(:,2),spsink{il},'natural');
ZZ = ZZ + f(XX,YY);
clear spos;
end
surf(XX,YY,ZZ)
shading interp
colormap(jet);
colorbar;
view(0,90)
axis off

hold on
alatL = rot_cell;
line([0 alatL(1,1)],[0 0],[1 1])
line([0 -alatL(2,1)],[0 -alatL(2,2)],[1 1])
line([-alatL(2,1) -alatL(2,1)+alatL(1,1)],[-alatL(2,2) -alatL(2,2)],[1 1])
line([-alatL(2,1)+alatL(1,1) alatL(1,1)],[-alatL(2,2) 0],[1 1])
C = [1/2*alatL(1,1),-1/3*alatL(2,2),1.1] ;   % center of circle
R = 1.5;
teta=0:0.01:2*pi ;
x=C(1)+R*cos(teta);
y=C(2)+R*sin(teta) ;
z = C(3)+zeros(size(x)) ;
patch(x,y,z,'g')
plot3(C(1),C(2),C(3),'k')
C = [0,-2/3*alatL(2,2),1.1] ;   % center of circle
R = 1.5;
teta=0:0.01:2*pi ;
x=C(1)+R*cos(teta);
y=C(2)+R*sin(teta) ;
z = C(3)+zeros(size(x)) ;
patch(x,y,z,'c')
plot3(C(1),C(2),C(3),'k')
C = [-alatL(2,1),-alatL(2,2),1.1] ;   % center of circle
R = 1.5;
teta=0:0.01:2*pi ;
x=C(1)+R*cos(teta);
y=C(2)+R*sin(teta) ;
z = C(3)+zeros(size(x)) ;
patch(x,y,z,'m')
plot3(C(1),C(2),C(3),'k')
%scatter3(supercell(:,1),supercell(:,2),supercell(:,3),10,spsink(:),'filled')


clear XX YY ZZ x y rot_cell supercell spsink
saveas(gca,join([filename,'_psink_',k,'.fig']))
end
