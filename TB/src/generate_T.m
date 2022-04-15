%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to generate rotation matrix from 
% symmetrised to unsymmetrised basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% orbitals - orbitals class
% norbs - Number of orbitals
%
% Outputs:
% T - rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = generate_T(orbitals,norbs,nspin,nlayer)
T = zeros(norbs);
for spin = -nspin + 1 : 2 : nspin - 1
   orb1 = findobj(orbitals,'l',2,'Spin',spin);
   diorb = [orb1.Ham_index];
   clear orb1;
   for ii = 1 : size(diorb,2)
       iorb = diorb(ii);
       for jj = 1 : size(diorb,2)
          jorb = diorb(jj);
          if((orbitals(jorb).Centre == orbitals(iorb).Centre))
             if (orbitals(jorb).Rel_index == orbitals(iorb).Rel_index)
                T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1.0;
             end
          end
       end
   end

   for ilayer = 1 : nlayer
   orb1 = findobj(orbitals,'l',1,'Spin', spin,'Layer',ilayer);
   piorb = [orb1.Ham_index];
   clear orb1;
   for ii = 1 : size(piorb,2)
      iorb = piorb(ii);
      for jj = 1 : size(piorb,2)
         jorb = piorb(jj);
         if(norm(orbitals(jorb).Centre-orbitals(iorb).Centre) < 0.1)
            if(orbitals(iorb).Rel_index == 3 && orbitals(jorb).Rel_index == 3)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 4 && orbitals(jorb).Rel_index == 4)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = -1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 5 && orbitals(jorb).Rel_index == 5)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = -1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 9 && orbitals(jorb).Rel_index == 9)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 10 && orbitals(jorb).Rel_index == 10)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 11 && orbitals(jorb).Rel_index == 11)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 3 && orbitals(jorb).Rel_index == 9)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 4 && orbitals(jorb).Rel_index == 10)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 5 && orbitals(jorb).Rel_index == 11)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 9 && orbitals(jorb).Rel_index == 3)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = -1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 10 && orbitals(jorb).Rel_index == 4)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            elseif(orbitals(iorb).Rel_index == 11 && orbitals(jorb).Rel_index == 5)
               T(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = 1/sqrt(2);
            end
         end
      end
   end
   end
end
        
end
