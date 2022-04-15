%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to find intralayer and interlayer 
% neighbors for each orbital
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% orbitals - orbitals class
% norbs - Number of orbitals
% cell - lattice vectors of moire cell (ang)
% ncell1,ncell2 - number of repeating cells in 
%                 a1 and a2 directions, respectively
% aXM -  minimum distance between a metal atom centre  
%        projected chalcogen atom centre onto metal plane
% nspin - npsins = 2 with SOC, nspin = 1 otherwise
%
% Outputs:
% orbitals - orbitals class 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale
% dz2-pz neighbor list added by Kemal Atalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orbitals] = opt_find_neigh2(orbitals,mcell,ncell1,ncell2,...
     bondlength,aXX2,nspin,type1,type2,type3,type4,type5,type6,...
     ia1,ia2,Rot,multilayer,nlayer,interlayer_int,Interpd,flipped)

% Make this an input later, X-M interlayer coupling cutoff                                                             
aXM = 6;

strain=0.1;
tol=strain*bondlength;

% Find nearest- 2nd-nearest- and 3rd-nearest- neighbors for each atom
for spin = -nspin + 1 : 2 : nspin - 1
    for layer = 1 : nlayer
        if(flipped(layer)==0)
            a1 = ia1*Rot{layer}';
            a2 = ia2*Rot{layer}';
        elseif(flipped(layer)==1)
            a1 = ia1*Rot{layer}';
            a2 = [ia2(1),-ia2(2),ia2(3)]*Rot{layer}';
        end
        for isym = [-1,1]
            orb1 = findobj(orbitals,'Layer',layer,'Spin',spin,'M1sym',isym);
            norbs = size(orb1,1);
            coor = [reshape([orb1(:).Centre],[3,norbs])'];
            hindex = [[orb1(:).Ham_index]'];
            supercell = zeros(ncell1^2*norbs,2);
            super_hindex = zeros(ncell1^2*norbs,3);
            dist = zeros(norbs,ncell1^2*norbs);
            ir = 0;
            for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                    ir = ir + 1;
                    supercell((ir-1)*norbs+1:ir*norbs,1:2) = coor(:,1:2) + n1*mcell(1,1:2)+n2*mcell(2,1:2);
                    super_hindex((ir-1)*norbs+1:ir*norbs,:) = [hindex(:),n1*ones(norbs,1),n2*ones(norbs,1)];
                end
            end
            dist = pdist2(coor(:,1:2),supercell);
            
            
            for ii = 1:norbs
                iorb = hindex(ii);
                % Find nearest neighbors
                dist_list = find(abs(dist(ii,:) - bondlength(layer)) <= tol(layer));
                orbitals(iorb).Intra_n_list = super_hindex(dist_list);
                if (orbitals(iorb).l == 2 && length(dist_list) > 9)
                    [size(dist_list,2),iorb]
                    error('d orbital: Too many neighbors')
                elseif (orbitals(iorb).l == 1 && length(dist_list) > 9)
                    [size(dist_list,2),iorb]
                    error('p orbital: Too many neighbors')
                end
                n1 = super_hindex(dist_list,2);
                n2 = super_hindex(dist_list,3);
                orbitals(iorb).Intra_n_centre = reshape([orbitals(super_hindex(dist_list)).Centre],...
                    [3,length(dist_list)])' + n1.*mcell(1,:) + n2.*mcell(2,:);
                delta = orbitals(iorb).Intra_n_centre - orbitals(iorb).Centre;
                M2sym = orbitals(iorb).M2sym*[orbitals(super_hindex(dist_list)).M2sym];
                ir = orbitals(iorb).Rel_index;
                jr = [orbitals(orbitals(iorb).Intra_n_list).Rel_index];
                %orbitals(iorb).Intra_n_hoppings = assign_hoppings(a1,a2,delta,1,layer,...
                %    M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
		loc_theta1 = orbitals(iorb).loc_theta;
		loc_theta2 = [orbitals(orbitals(iorb).Intra_n_list).loc_theta];
                orbitals(iorb).Intra_n_hoppings = assign_hoppings_all(a1,a2,delta,loc_theta1,...
		    loc_theta2,1,layer,M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
                
                % Find next-nearest neighbors
                dist_list = find(abs(dist(ii,:) - sqrt(3)*bondlength(layer)) <= sqrt(3)*tol(layer)); % *sqrt(3)) < sqrt(3)*tol(layer));
                orbitals(iorb).Intra_nn_list = super_hindex(dist_list);
                if (orbitals(iorb).l == 2 && length(dist_list) > 18)
                    error('d orbital: Too many neighbors')
                elseif (orbitals(iorb).l == 1 && length(dist_list) > 18)
                    error('p orbital: Too many neighbors')
                end
                n1 = super_hindex(dist_list,2);
                n2 = super_hindex(dist_list,3);
                orbitals(iorb).Intra_nn_centre = reshape([orbitals(super_hindex(dist_list)).Centre],...
                    [3,length(dist_list)])' + n1.*mcell(1,:) + n2.*mcell(2,:);
                delta = orbitals(iorb).Intra_nn_centre - orbitals(iorb).Centre;
                M2sym = orbitals(iorb).M2sym*[orbitals(super_hindex(dist_list)).M2sym];
                ir = orbitals(iorb).Rel_index;
                jr = [orbitals(orbitals(iorb).Intra_nn_list).Rel_index];
                %orbitals(iorb).Intra_nn_hoppings = assign_hoppings(a1,a2,delta,2,layer,...
                %    M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
		loc_theta1 = orbitals(iorb).loc_theta;
		loc_theta2 = [orbitals(orbitals(iorb).Intra_nn_list).loc_theta];
                orbitals(iorb).Intra_nn_hoppings = assign_hoppings_all(a1,a2,delta,loc_theta1,...
		    loc_theta2,2,layer,M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
                
                % Find next-next-nearest neighbors
                dist_list = find(abs(dist(ii,:) - 2.0*bondlength(layer)) <= 2.0*tol(layer)); % < 2.0*tol(layer));
                orbitals(iorb).Intra_nnn_list = super_hindex(dist_list);
                if (orbitals(iorb).l == 2 && length(dist_list) > 9)
                    error('d orbital: Too many neighbors')
                elseif (orbitals(iorb).l == 1 && length(dist_list) > 9)
                    error('p orbital: Too many neighbors')
                end
                n1 = super_hindex(dist_list,2);
                n2 = super_hindex(dist_list,3);
                orbitals(iorb).Intra_nnn_centre = reshape([orbitals(super_hindex(dist_list)).Centre],...
                    [3,length(dist_list)])' + n1.*mcell(1,:) + n2.*mcell(2,:);
                if(isym > 0)
                    delta = orbitals(iorb).Intra_nnn_centre - orbitals(iorb).Centre;
                    M2sym = orbitals(iorb).M2sym*[orbitals(super_hindex(dist_list)).M2sym];
                    ir = orbitals(iorb).Rel_index;
                    isym1_jr = findobj(orbitals(super_hindex(dist_list)),'M1sym',isym);
                    jr = [isym1_jr.Rel_index];
                    %orbitals(iorb).Intra_nnn_hoppings = assign_hoppings(a1,a2,delta,3,layer,...
                    %    M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
	            loc_theta1 = orbitals(iorb).loc_theta;
		    loc_theta2 = [isym1_jr.loc_theta];
                    orbitals(iorb).Intra_nnn_hoppings = assign_hoppings_all(a1,a2,delta,loc_theta1,...
	                loc_theta2,3,layer,M2sym,ir,jr,type1,type2,type3,type4,type5,type6);
                    clear M2sym ir jr delta
                end
            end
        end
    end
    clear dist dist_list coor orb1 supercell norbs
    if(multilayer && interlayer_int)
        for layer = 1 : nlayer -1
            orb1 = findobj(orbitals,'Layer',layer,'Spin',spin,'l',1);
            orb2 = findobj(orbitals,'Layer',layer+1,'Spin',spin,'l',1);
            norbs1 = size(orb1,1);
            norbs2 = size(orb2,1);
            coor1 = [reshape([orb1(:).Pcentre],[3,norbs1])'];
            coor2 = [reshape([orb2(:).Pcentre],[3,norbs2])'];
            hindex1 = [[orb1(:).Ham_index]'];
            hindex2 = [[orb2(:).Ham_index]'];
            ir = 0;
            supercell = zeros(ncell1^2*norbs2,3);
            super_hindex2 = zeros(ncell1^2*norbs2,3);
            for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                    ir = ir + 1;
                    supercell((ir-1)*norbs2+1:ir*norbs2,:) = coor2 + n1*mcell(1,:)+n2*mcell(2,:);
                    super_hindex2((ir-1)*norbs2+1:ir*norbs2,:) = [hindex2(:),n1*ones(norbs2,1),n2*ones(norbs2,1)];
                end
            end
            dist = pdist2(coor1,supercell);
            for ii = 1 : norbs1
                iorb = hindex1(ii);
                % Find nearest neighbor
                dist_list = find(dist(ii,:) <= aXX2);
                n1 = super_hindex2(dist_list,2);
                n2 = super_hindex2(dist_list,3);
                orbitals(iorb).Inter_nn_list = super_hindex2(dist_list);
                orbitals(iorb).Inter_nn_centre = reshape([orbitals(super_hindex2(dist_list)).Pcentre],[3,length(dist_list)])' - orbitals(iorb).Pcentre + n1.*mcell(1,:) + n2.*mcell(2,:);
            end
            clear orb1 orb2
            % Interlayer X-M neighbours
            orb1 = findobj(orbitals,'Layer',layer,'Spin',spin,{'Rel_index',9,'-or','Rel_index',6,'-or','Rel_index',7,'-or','Rel_index',8});
            orb2 = findobj(orbitals,'Layer',layer+1,'Spin',spin,{'Rel_index',3,'-or','Rel_index',6,'-or','Rel_index',7,'-or','Rel_index',8});
            diorb = [orb1(:).Ham_index];
            djorb = [orb2(:).Ham_index];
            clear orb2;
            clear orb1;
            s1 = size(diorb,2);
            s2 = size(djorb,2);
            for ii = 1 : s1
                iorb = diorb(ii);
                in4 = 0;
                in4_list = zeros(1);
                in4_centre = zeros(1,3);
                for jj = 1 : s2
                    jorb = djorb(jj);
                    if((orbitals(iorb).Rel_index == 9 && orbitals(jorb).Rel_index ~= 3) || (orbitals(iorb).Rel_index ~= 9 && orbitals(jorb).Rel_index == 3))
                        % Loop over all sites in the supercell
                        for n1 = -(ncell1-1)/2 : (ncell1-1)/2
                            for n2 = -(ncell2-1)/2 : (ncell2-1)/2
                                if(orbitals(iorb).Rel_index == 9)
                                    shift = n1*mcell(1,:)+n2*mcell(2,:) + orbitals(jorb).Centre;
                                    d = shift - orbitals(iorb).Pcentre;
                                else
                                    shift = n1*mcell(1,:)+n2*mcell(2,:) + orbitals(jorb).Pcentre;
                                    d = shift - orbitals(iorb).Centre;
                                end
                                dsqr = norm(d);
                                % Find neighbors in given cutoff radius
                                if(dsqr < aXM)
                                    in4 = in4 + 1;
                                    in4_list(in4) = orbitals(jorb).Ham_index;
                                    in4_centre(in4,:) = d;
                                end
                            end
                        end
                    end
                end
                orbitals(iorb).Inter_mx_n_list = in4_list;
                orbitals(iorb).Inter_mx_n_centre = in4_centre;
            end
        end
    end
end

clear in_list in_centre inn_list inn_centre innn_list innn_centre in2_list in2_centre in3_list in3_centre dis1 dis2
clear liorb ljorb pdiorb pdjorb dpiorb dpjorb diorb djorb  
end
