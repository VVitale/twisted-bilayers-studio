%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to add SOC to TB Hamiltonian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% Hmat - Hamiltonian matrix
% Hsoc - 22x22 SOC matrix in 22-orbital basis
% orbitals - orbitals class
% norbs - Number of orbitals
%
% Outputs:
% Hmat - Hamiltonian matrix with SOC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Hmat = add_soc(Hmat,Hsoc,orbitals,ilayer)

liorb = cell(2,2);

for ispin = 1 : 2
    if(ispin == 1)
        iispin = -1;
    else
        iispin = 1;
    end
    for il = 1 : 2
         orb = findobj(orbitals,'Spin',iispin,'l',il,'Layer',ilayer);
         liorb{ispin,il} = [orb.Ham_index];
    end
end
clear orb
orb1 = liorb{1,2};
orb2 = liorb{2,2};
orb3 = liorb{1,1};
orb4 = liorb{2,1};
for ii = 1 : length(orb1)
    iorb = orb1(ii);
    for jj = 1 : length(orb1)
        jorb = orb1(jj);
        if(orbitals(jorb).Centre == orbitals(iorb).Centre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index);
        end
    end
    for jj = 1 : length(orb2)
        jorb = orb2(jj);
        if(orbitals(jorb).Centre == orbitals(iorb).Centre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index+11);
            Hmat(orbitals(jorb).Ham_index,orbitals(iorb).Ham_index) = Hmat(orbitals(jorb).Ham_index,orbitals(iorb).Ham_index) + conj(Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index+11));% + conj(Hsoc(orb1.Rel_index,orb2.Rel_index));
        end
    end
end
for ii = 1 : length(orb2)
    iorb = orb2(ii);
    for jj = 1 : length(orb2)
        jorb = orb2(jj);
        if(orbitals(jorb).Centre == orbitals(iorb).Centre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index+11,orbitals(jorb).Rel_index+11);
        end
    end
end

for ii = 1 : length(orb3)
    iorb = orb3(ii);
    for jj = 1 : length(orb3)
        jorb = orb3(jj);
        if(orbitals(jorb).Pcentre == orbitals(iorb).Pcentre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index);
        end
    end
    for jj = 1 : length(orb4)
        jorb = orb4(jj);
        if(orbitals(jorb).Pcentre == orbitals(iorb).Pcentre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index+11);
            Hmat(orbitals(jorb).Ham_index,orbitals(iorb).Ham_index) = Hmat(orbitals(jorb).Ham_index,orbitals(iorb).Ham_index) + conj(Hsoc(orbitals(iorb).Rel_index,orbitals(jorb).Rel_index+11));% + conj(Hsoc(orb1.Rel_index,orb2.Rel_index));
            %Hmat(orbitals(jorb).Ham_index,orbitals(iorb).Ham_index) = conj(Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index));% + conj(Hsoc(orb1.Rel_index,orb2.Rel_index));
        end
    end
end


for ii = 1 : length(orb4)
    iorb = orb4(ii);
    for jj = 1 : length(orb4)
        jorb = orb4(jj);
        if(orbitals(jorb).Pcentre == orbitals(iorb).Pcentre)
            Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) = Hmat(orbitals(iorb).Ham_index,orbitals(jorb).Ham_index) + Hsoc(orbitals(iorb).Rel_index+11,orbitals(jorb).Rel_index+11);
        end
    end
end
clear orb1 orb2 orb3 orb4 liorb
end
