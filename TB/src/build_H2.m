%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to compute Hamiltonian matrix
% This routine generalises the TB Hamiltonian
% in PRB 92, 205108 (2015) for different systems
% This routines generates the TB Hamiltonian
% matrix for monolayers and bilayers. 
% In the case of bilayers this routines computes 
% the block-diagonal part for each layer.
% The off-diagonal interlayer interaction is 
% computed in a different routine Vint.m
% In the case of a SOC calculation this routines 
% only compute the block-diagonal part of the 
% Hamiltonian for spin-up and spin-down channel
% The off-diagonal blocks are added by add_soc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% Hmat  - Hamiltonian matrix
% norbs - Number of orbitals
% orbitals - orbitals class
% ia1,ia2,ia3 - Planar vectors used to assign hoppings
% k - k-point in Cartesian coordinates
% onsite,type1,type2,type3,type4,type5,type6 - All TB parameters for monolayer
% Rot - Rotation matrix for upper layer
% bilayer - Whether dealing with a bilayer or monolayer
% nspins - nspins = 2 for SOC, otherwise nspins = 1
% interlayer_int - whether to compute interlayer interaction
% lsoc - Whether to add SOC
% T - unsymmetrising matrix
% Vsigma - sigma parameter for p-orbitals in
%          Slater-Koster rule
% Vpi    - pi parameter for p-orbitals in 
%          Slater-Koster rule
% R,eta  - Parameters for spatial-dependent part 
%          in Slater-Koster rule [1]
% theta  - twist angle
% Hsoc   - SOC matrices for p and d orbitals
%
% Outputs:
% Hmat - Hamiltonian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hmat,gradH_x] = build_H2(Hmat,orbitals,Mxz,k,nlayer,multilayer,...
                        nspins,interlayer_int,lsoc,T,pp_vint_z,pp_vint_parm, ...
			pd_vint_lay1_z,pd_vint_lay1_parm,pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,...
                        Interpd,gradient,varargin)

if(gradient && isempty(varargin))
   error('Error in build_H.m: Hamiltonian gradient required. Aborting ...')
end

if(~gradient && ~isempty(varargin))
   error('Error in build_H.m: Found gradient but gradient not required. Aborting ...')
end

if(gradient && length(varargin)==1)
   gradH_x = cell2mat(varargin(1));
end

for layer = 1 : nlayer
    orb = findobj(orbitals,'Layer',layer);
    liorb = [orb.Ham_index];
    clear orb1;
    for ii = 1 : size(liorb,2)
        iorb = liorb(ii);
        if(orbitals(iorb).M1sym > 0)
            nnn = size(orbitals(iorb).Intra_nnn_list,2);
            % Contribution from next nearest neighbors
            for in = 1 : nnn
                jorb = orbitals(iorb).Intra_nnn_list(in);
                hopping = orbitals(iorb).Intra_nnn_hoppings(in);
                delta = orbitals(iorb).Intra_nnn_centre(in,1:3) - orbitals(iorb).Centre;
                
                % Find type of hopping
                % Take only even states
                hij = hopping*exp(1i*dot(delta,k));
                Hmat(iorb,jorb) = Hmat(iorb,jorb) + hij; 
                Hmat(jorb,iorb) = Hmat(jorb,iorb) + conj(hij); %hopping*exp(1i*dot(delta,k));
                if(gradient)
                    gradH_x(iorb,jorb) = gradH_x(iorb,jorb) - delta(1)*hij - 1i*delta(2)*hij;
                    gradH_x(jorb,iorb) = gradH_x(jorb,iorb) + conj(- delta(1)*hij - 1i*delta(2)*hij); %hopping*exp(1i*dot(delta,k));
                end
            end
        end
        nn1 = size(orbitals(iorb).Intra_nn_list,2);
        for in = 1 : nn1
            jorb = orbitals(iorb).Intra_nn_list(in);
            if(orbitals(jorb).Rel_index > orbitals(iorb).Rel_index)
                hopping = orbitals(iorb).Intra_nn_hoppings(in);
                delta = orbitals(iorb).Intra_nn_centre(in,1:3) - orbitals(iorb).Centre;            
                hij = hopping*exp(1i*dot(delta,k));
                Hmat(iorb,jorb) = Hmat(iorb,jorb) + hij;
                Hmat(jorb,iorb) = Hmat(jorb,iorb) + conj(hij); %hopping*exp(1i*dot(delta,k));
                if(gradient)
                   gradH_x(iorb,jorb) = gradH_x(iorb,jorb) - delta(1)*hij - 1i*delta(2)*hij;
                   gradH_x(jorb,iorb) = gradH_x(jorb,iorb) + conj(- delta(1)*hij - 1i*delta(2)*hij); %hopping*exp(1i*dot(delta,k));
                end
            end
        end
        n = size(orbitals(iorb).Intra_n_list,2);
        % Contribution from nearest neighbors
        for in = 1 : n
            jorb = orbitals(iorb).Intra_n_list(in);
            hopping = orbitals(iorb).Intra_n_hoppings(in);
            delta = orbitals(iorb).Intra_n_centre(in,1:3) - orbitals(iorb).Centre;
                        
            % Find type of hopping
            hij = hopping*exp(1i*dot(delta,k));
            Hmat(iorb,jorb) = Hmat(iorb,jorb) + hij;
            if(gradient)
                gradH_x(iorb,jorb) = gradH_x(iorb,jorb) - delta(1)*hij - 1i*delta(2)*hij;
            end
        end
        % Contribution from next nearest neighbors
        nn2 = size(orbitals(iorb).Intra_nn_list,2);
        for in = 1 : nn2
            jorb = orbitals(iorb).Intra_nn_list(in);
            if(orbitals(jorb).Rel_index == orbitals(iorb).Rel_index)
                hopping = orbitals(iorb).Intra_nn_hoppings(in);
                delta = orbitals(iorb).Intra_nn_centre(in,1:3) - orbitals(iorb).Centre;
                
                % Get connecting vector between centres
                
                % Find type of hopping
                hij = hopping*exp(1i*dot(delta,k));
                Hmat(iorb,jorb) = Hmat(iorb,jorb) + hij;
                if(gradient)
                    gradH_x(iorb,jorb) = gradH_x(iorb,jorb) - delta(1)*hij - 1i*delta(2)*hij;
                end
            end
        end
    end
end

H_rotated = false;
% Add interlayer interaction
if(multilayer && interlayer_int)
    % Rotate Hamiltonian into unsymmetrised basis
    Hmat = T'*(Hmat*T);
    % Apply mask for 2H stackings
    Hmat = Mxz'*Hmat*Mxz;
    % Check Hamiltonian is Hermitian
    if(gradient)
       gradH_x = T'*(gradH_x*T);
       [Hmat,gradH_x] = Vint(Hmat,nlayer,orbitals,pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm,pd_vint_lay2_z,pd_vint_lay2_parm,k,theta,Interpd,nspins,gradient,gradH_x);
    else
       [Hmat] = Vint(Hmat,nlayer,orbitals,pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm,pd_vint_lay2_z,pd_vint_lay2_parm,k,theta,Interpd,nspins,gradient);
    end
    H_rotated = true;
end

% Add SOC
if (lsoc)
    if (~H_rotated)
       Hmat = T'*(Hmat*T);
       H_rotated = true;
    end
    for ilayer = 1 : nlayer
       Hmat = add_soc(Hmat,Hsoc{ilayer},orbitals,ilayer);
    end
end


%check if Hamiltonian is Hermitian
%if(normest(Hmat-Hmat')>1E-4)
%   error('Hamiltonian is not hermitian!')
%end
Hmat = 0.5*(Hmat+Hmat');
if(gradient)
   gradH_x = 0.5*(gradH_x+gradH_x');
end
end
