%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine sets properties for of each      %
% orbital, e.g. Hamiltonian index, Relative     %
% index, symmetry eigenvalue wrt to M1 and M2   %
% etc.                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orbitals] = initialise_orbitals(orbitals,structure,mcell,tot_natoms,natoms,multilayer,nlayer,nspin)
jorb = 0;
layer1 = 1;
if (multilayer)
    layer2 = nlayer;
else
    layer2 = layer1;
end

% VV: TODO
% Principal quantum numbers and atomic numbers
% for ilayer = 1 : nlayer
% if(any(structure.name(1:natoms(ilayer)) == "Mo"))
%     metal_pqm(ilayer) = 4;
%     ZM(ilayer) = 42; 
% elseif(any(structure.name(1:natoms(ilayer)) == "W"))
%     metal_pqm(ilayer) = 5;
%     ZM(ilayer) = 74;
% end
% 
% if(any(structure.name(1:natoms(ilayer)) == "S"))
%     chalc_pqm(1) = 3;
%     ZX(1) = 16;
% elseif(any(structure.name(1:natoms(1)) == "Se"))
%     chalc_pqm(1) = 4;
%     ZX(1) = 34;
% end
% 
% if(any(structure.name(natoms(1)+1:tot_natoms) == "Mo"))
%     metal_pqm(2) = 4;
%     ZM(2) = 42; 
% elseif(any(structure.name(natoms(1)+1:tot_natoms) == "W"))
%     metal_pqm(2) = 5;
%     ZM(2) = 74;
% end
% 
% if(any(structure.name(natoms(1)+1:tot_natoms) == "S"))
%     chalc_pqm(2) = 3;
%     ZX(2) = 16;
% elseif(any(structure.name(natoms(1)+1:tot_natoms) == "Se"))
%     chalc_pqm(2) = 4;
%     ZX(2) = 34;
% end
% end


isfree(1:tot_natoms) = true;
finish = 0;
% Loop over atoms
for layer = layer1 : layer2
    start = finish + 1;
    finish = start + natoms(layer) -1;
    for spin = -nspin+1 : 2 : nspin - 1
        isfree(start:finish) = true;
        for inat = start : finish
            % Find a free metal atom
            if ((structure.name(inat) == "Mo" || structure.name(inat) == "W") && structure.layer(inat) == layer && isfree(inat))
                isfree(inat) = false;
                x1 = [structure.x(inat),structure.y(inat),structure.z(inat)];
                for iorb = 1:11
                    if (iorb == 1 || iorb == 2 || iorb == 6 || iorb == 7 || iorb == 8)
                        jorb = jorb + 1;
                        orbitals(jorb).Ham_index = jorb;
                        orbitals(jorb).Rel_index = iorb;
                        % Associate centre to orbitals
                        % Associate a Metal atom to this orbital
                        orbitals(jorb).Found_centre = true;
                        orbitals(jorb).Centre = x1;
                        orbitals(jorb).Layer = structure.layer(inat);
                        orbitals(jorb).Spin = spin;
                        %orbitals(jorb).N = metal_pqm(layer);
                        %orbitals(jorb).Z = ZM(layer);
                    end
                    
                end
            end
        end
        
        for inat = start : finish
            foundX = false;
            if ((structure.name(inat) == "S" || structure.name(inat) == "Se") && structure.layer(inat) == layer && isfree(inat))
                incell = 0;
                x1 = [structure.x(inat),structure.y(inat),structure.z(inat)];
                isfree(inat) = false;
                % find its companion
                for jnat = start : finish
                    if((structure.name(jnat) == "S" || structure.name(jnat) == "Se") && structure.layer(jnat) == layer && isfree(jnat) )
                        
                        x2 = [structure.x(jnat),structure.y(jnat),structure.z(jnat)];
                        for in = -1 : 1
                            for im = -1 : 1
                                shift = x2 + in*mcell(1,:) + im*mcell(2,:);
                                dsqr = norm(x1(1:2) - shift(1:2));
                                if (dsqr<0.5)
                                    foundX = true;
                                    isfree(jnat) = false;
                                    incell = incell + 1;
                                    if (incell > 1)
                                        error('More than 2 atoms per unit cell')
                                    elseif (incell == 0)
                                        error('No neighbors in the same cell found')
                                    end
                                    break
                                end
                            end
                        end
                        if(foundX)
                            break;
                        end
                    end
                end
                
                if (foundX)
                    for iorb = 1 : 11
                        if (iorb == 3 || iorb == 4 || iorb == 5)
                            jorb = jorb + 1;
                            orbitals(jorb).Ham_index = jorb;
                            orbitals(jorb).Rel_index = iorb;
                            orbitals(jorb).Found_centre = true;
                            if(x1(3)<x2(3))
                                orbitals(jorb).Centre = [x1(1),x1(2),(x1(3)+x2(3))/2];
                                orbitals(jorb).Pcentre = [x1(1),x1(2),x1(3)];
                            else
                                orbitals(jorb).Centre = [x2(1),x2(2),(x1(3)+x2(3))/2];
                                orbitals(jorb).Pcentre = [x2(1),x2(2),x2(3)];
                            end
                            orbitals(jorb).Layer = structure.layer(inat);
                            orbitals(jorb).Spin = spin;
                            %orbitals(jorb).N = chalc_pqm(layer);
                            %orbitals(jorb).Z = ZX(layer);
%                             if(layer~=1 && layer~=nlayer)
%                                 orbitals(jorb).Inner = true;
%                             elseif(layer==1)
%                                orbitals(jorb).Inner = false;
%                             elseif(layer==nlayer)
%                                 orbitals(jorb).Inner = true;
%                             end
                            orbitals(jorb).Inner = false;
                        elseif(iorb == 9 || iorb == 10 || iorb == 11)
                            jorb = jorb + 1;
                            orbitals(jorb).Ham_index = jorb;
                            orbitals(jorb).Rel_index = iorb;
                            orbitals(jorb).Found_centre = true;
                            if(x1(3)>x2(3))
                                orbitals(jorb).Centre = [x1(1),x1(2),(x1(3)+x2(3))/2];
                                orbitals(jorb).Pcentre = [x1(1),x1(2),x1(3)];
                            else
                                orbitals(jorb).Centre = [x2(1),x2(2),(x1(3)+x2(3))/2];
                                orbitals(jorb).Pcentre = [x2(1),x2(2),x2(3)];
                            end
                            orbitals(jorb).Layer = structure.layer(inat);
                            orbitals(jorb).Spin = spin;
                            %orbitals(jorb).N = chalc_pqm(layer);
                            %orbitals(jorb).Z = ZX(layer);
%                             if(layer~=1 && layer~=nlayer)
%                                 orbitals(jorb).Inner = true;
%                             elseif(layer==1)
%                                 orbitals(jorb).Inner = true;
%                             elseif(layer==nlayer)
%                                 orbitals(jorb).Inner = false;
%                             end
                            if(layer~=nlayer)
                                orbitals(jorb).Inner = true;
                            else
                                orbitals(jorb).Inner = false;
                            end
                        end
                    end
                elseif(~foundX && jnat~=natoms(layer) && jnat~=tot_natoms)
                    jnat
                    error('Error in initialise_orbitals.m: Could not find two chalcogens in the same cell')
                end
            end
        end
    end
end
list = [orbitals.Ham_index];
ind = find(list==0);
if(~isempty(ind))
   ind
   error('Error in initialise_orbitals.m: Cannot initialise orbitals')
end
clear orb1 list
end
