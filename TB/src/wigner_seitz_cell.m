function [structure_in_ws, id_in_ws, IN,tr_vec] = wigner_seitz_cell(convention,structure,mcell,nat)

% Generate 3x3 supercell
supercell = zeros(4*nat,3);
superid = zeros(4*nat,1);
structure_in_ws = zeros(size(supercell));
id_in_ws = zeros(size(superid));
is = 0;
min_dis = 100;
for im = [0,-1]
   for jm = [0,-1]
      for inat = 1 : nat
         is = is + 1;
         supercell(is,:) = [structure.x(inat),structure.y(inat),structure.z(inat)] + im*mcell(1,:) + jm*mcell(2,:);
	 if(norm(supercell(is,:)) < min_dis && (strcmp(structure.name(inat),'Mo') || strcmp(structure.name(inat),'W')))
		 min_dis = norm(supercell(is,:));
		 min_at_id = is;
	 end
         superid(is) = inat;
      end
   end
end
% Translate supercell such that a metal atom is at the origin.
tr_vec = supercell(min_at_id,:);
supercell = supercell - tr_vec;

% Define WS cell
delta = 0.0;
if(convention==2)
   % Define vertices of polyong
   x1 = 1/3*(mcell(1,1:2) + mcell(2,1:2));
   x2 = 1/3*(-mcell(1,1:2) + 2*mcell(2,1:2));
   x3 = 1/3*(-2*mcell(1,1:2) + mcell(2,1:2));
   x4 =-1/3*(mcell(1,1:2) + mcell(2,1:2));
   x5 = 1/3*(mcell(1,1:2) - 2*mcell(2,1:2));
   x6 = 1/3*(2*mcell(1,1:2) - mcell(2,1:2));
elseif(convention==1)
   % Define vertices of polyong
   x1 = 1/3*(2*mcell(1,1:2) + mcell(2,1:2));
   x2 = 1/3*(mcell(1,1:2) + 2*mcell(2,1:2));
   x3 = 1/3*(-mcell(1,1:2) + mcell(2,1:2));
   x4 =-1/3*(2*mcell(1,1:2) + mcell(2,1:2));
   x5 =-1/3*(mcell(1,1:2) + 2*mcell(2,1:2));
   x6 = 1/3*(mcell(1,1:2) - mcell(2,1:2));
end
XV = [x1(1);x2(1);x3(1);x4(1)-delta;x5(1)-delta;x6(1)-delta;x1(1)];
YV = [x1(2);x2(2);x3(2);x4(2)-delta;x5(2)-delta;x6(2)-delta;x1(2)];
IN = inpolygon(supercell(:,1),supercell(:,2),XV,YV);
structure_in_ws = supercell(IN,:);
id_in_ws = superid(IN);
% Check number of atoms in WS cell is the same as moire cell
nat_ws = size(structure_in_ws,1);
if(nat_ws ~= nat)
   nat_ws
   nat
   error('Number of atoms in WS cell is not the equal to Nat')
end

end
