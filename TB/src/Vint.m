%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to compute the off-diagonal 
% interlayer interaction of TB Hamiltonian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% Hmat - Hamiltonian matrix
% orbitals - orbitals class
% norbs - Number of orbitals
% Vsigma - sigma parameter for p-orbitals in
%          Slater-Koster rule
% Vpi    - pi parameter for p-orbitals in 
%          Slater-Koster rule
% R,eta  - Parameters for spatial-dependent part 
%          in Slater-Koster rule [1]
% k      - k-point in Cartesian coordinates
% theta  - twist angle
%
% Outputs:
% Hmat - Hamiltonian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] PRB 92, 205108 (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
% dz2-pz interaction added by Kemal Atalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hmat,gradH_x] = Vint(Hmat,nlayer,orbitals,pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
			  pd_vint_lay2_z,pd_vint_lay2_parm,k,theta,Interpd,nspin,gradient,varargin)
% Make this input later, fit to SK parameters for p-d hopping
checkhop=false;

if(gradient && isempty(varargin))
   error('Error in Vint.m: Hamiltonian gradient required. Aborting ...')
end

if(~gradient && ~isempty(varargin))
   error('Error in Vint.m: Found gradient in Vint.m but gradient not required. Aborting ...')
end 
if(gradient && length(varargin)==1)
   gradH_x = cell2mat(varargin(1));
end
% List to store hopping at distances - for checks of calculation
pdr = [];
pdt = [];

%   ct=cos(pi-theta);
%   st=sin(pi-theta);
for spin = -nspin+1:2:nspin-1
    for layer = 1 : nlayer-1
        orb1 = findobj(orbitals,'Layer',layer,'l',1,'Spin',spin,'Inner',true);
        liorb = [orb1(:).Ham_index];
        clear orb1;
        for ii = 1 : size(liorb,2)
            iorb = liorb(ii);
            if (size(orbitals(iorb).Inter_nn_list,2)>1)
                orb1 = orbitals(iorb);
                in = size(orb1.Inter_nn_list,2);
                for inn = 1 : in
                    orb2 = orbitals(orb1.Inter_nn_list(inn));
                    delta = orb1.Inter_nn_centre(inn,:);
                    dsqr = norm(delta);
                    nco = delta(3)/dsqr; mco = delta(2)/dsqr; lco = delta(1)/dsqr;
                    [Vpsig,Vppi] = interpolate_SK(1,dsqr,delta(3),pp_vint_z(:,layer),pp_vint_parm(:,layer));
                    if(orb1.Rel_index == 9 && orb2.Rel_index == 3)
                        hopping = set_hopping('pzpz',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 10 && orb2.Rel_index == 4)
                        hopping = set_hopping('pxpx',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 11 && orb2.Rel_index == 5)
                        hopping = set_hopping('pypy',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 10 && orb2.Rel_index == 3)
                        hopping = set_hopping('pxpz',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 11 && orb2.Rel_index == 3)
                        hopping = set_hopping('pypz',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 10 && orb2.Rel_index == 5)
                        hopping = set_hopping('pxpy',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 11  && orb2.Rel_index == 4)
                        hopping = set_hopping('pypx',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 9 && orb2.Rel_index == 5)
                        hopping = set_hopping('pzpy',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    elseif(orb1.Rel_index == 9  && orb2.Rel_index == 4)
                        hopping = set_hopping('pzpx',Vpsig,Vppi,lco,mco,nco,theta,layer);
                    else
                        error('No interlayer set')
                    end
                    Hmat(orb1.Ham_index,orb2.Ham_index) = Hmat(orb1.Ham_index,orb2.Ham_index) + hopping*exp(1i*dot(delta,k));
                    Hmat(orb2.Ham_index,orb1.Ham_index) = Hmat(orb2.Ham_index,orb1.Ham_index) + conj(hopping*exp(1i*dot(delta,k))); % + hopping*exp(1i*dot(delta,k));
                    if(gradient)
                        gradH_x(orb1.Ham_index,orb2.Ham_index) = gradH_x(orb1.Ham_index,orb2.Ham_index) + 1i*delta(1)*hopping*exp(1i*dot(delta,k));%- 1i*1i*delta(2)*hopping*exp(1i*dot(delta,k));
                        gradH_x(orb2.Ham_index,orb1.Ham_index) = gradH_x(orb2.Ham_index,orb1.Ham_index) + conj(1i*delta(1)*hopping*exp(1i*dot(delta,k)));%- 1i*1i*delta(2)*hopping*exp(1i*dot(delta,k)));% + 1i*delta(1)*hopping*exp(1i*dot(delta,k));
                    end
                end
            end
        end
        clear in;
        if (Interpd)
            orb1 = findobj(orbitals,'Spin',spin,'Layer',layer,{'Rel_index',9,'-or','Rel_index',6,'-or','Rel_index',7,'-or','Rel_index',8});
            liorb = [orb1(:).Ham_index];
            clear orb1;
            for ii = 1 : size(liorb,2)
                iorb = liorb(ii);
                if (size(orbitals(iorb).Inter_mx_n_list,2)>1)
                    orb1 = orbitals(iorb);
                    in = size(orb1.Inter_mx_n_list,2);
                    for inn = 1 : in
                        hopping = 0;
                        orb2 = orbitals(orb1.Inter_mx_n_list(inn));
                        delta = orb1.Inter_mx_n_centre(inn,:);
                        dsqr = norm(delta);
                        nco = delta(3)/dsqr; mco = delta(2)/dsqr; lco = delta(1)/dsqr;
                        if(orb1.Rel_index == 9 && orb2.Rel_index == 6  || orb1.Rel_index == 6 && orb2.Rel_index == 3)
                            % Determine the correct SK parameters
                            if(orb1.Rel_index == 6)
                                [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay1_z(:,layer),pd_vint_lay1_parm(:,layer));
                            else
                                [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay2_z(:,layer),pd_vint_lay2_parm(:,layer));
                            end
                            hopping = set_hopping('dz2pz',Vpdsig,Vpdpi,lco,mco,nco,theta,layer);
                            if(~ismembertol(dsqr,pdr,0.000001))
                                pdr = [pdr, dsqr];
                                pdt = [pdt, hopping];
                            end
                       elseif(orb2.Rel_index == 7 && orb1.Rel_index == 9)
                           % Determine the correct SK parameters
                           [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay2_z(:,layer),pd_vint_lay2_parm(:,layer));
                           % hopping = nco*(nco^2 - 0.5*(lco^2 + mco^2))*Vpdsig + sqrt(3)*nco*(lco^2 + mco^2)*Vpdpi;
                           hopping = set_hopping('pzdxy',Vpdsig,Vpdpi,lco,mco,nco,theta,layer);
                       elseif(orb1.Rel_index == 7 && orb2.Rel_index == 3)
                           [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay1_z(:,layer),pd_vint_lay1_parm(:,layer));
                           hopping = set_hopping('dxypz',Vpdsig,Vpdpi,lco,mco,nco,theta,layer);
                       elseif(orb2.Rel_index == 8 && orb1.Rel_index == 9)
                           % Determine the correct SK parameters
                           [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay2_z(:,layer),pd_vint_lay2_parm(:,layer));
                           hopping = set_hopping('pzdx2y2',Vpdsig,Vpdpi,lco,mco,nco,theta,layer);
                           %hopping = sqrt(3)*nco*lco*mco*Vpdsig - 2*nco*lco*mco*Vpdpi;
                       elseif(orb1.Rel_index == 8 && orb2.Rel_index == 3)
                          [Vpdsig,Vpdpi] = interpolate_SK(2,dsqr,delta(3),pd_vint_lay1_z(:,layer),pd_vint_lay1_parm(:,layer));
                          hopping = set_hopping('dx2y2pz',Vpdsig,Vpdpi,lco,mco,nco,theta,layer);
                        end
                        
                        Hmat(orb1.Ham_index,orb2.Ham_index) = Hmat(orb1.Ham_index,orb2.Ham_index) + hopping*exp(1i*dot(delta,k));
                        Hmat(orb2.Ham_index,orb1.Ham_index) = Hmat(orb2.Ham_index,orb1.Ham_index) + conj(hopping*exp(1i*dot(delta,k)));% + hopping*exp(1i*dot(delta,k));
                        if(gradient)
                            gradH_x(orb1.Ham_index,orb2.Ham_index) = gradH_x(orb1.Ham_index,orb2.Ham_index) + 1i*delta(1)*hopping*exp(1i*dot(delta,k));%- 1i*1i*delta(2)*hopping*exp(1i*dot(delta,k));
                            gradH_x(orb2.Ham_index,orb1.Ham_index) = gradH_x(orb2.Ham_index,orb1.Ham_index) + conj(1i*delta(1)*hopping*exp(1i*dot(delta,k)));%- 1i*1i*delta(2)*hopping*exp(1i*dot(delta,k)));% + 1i*delta(1)*hopping*exp(1i*dot(delta,k));
                        end
                    end
                end
            end
        end
    end
end
if(checkhop && norm(k) == 0)
    pdr
    pdt
end
clear orb1 orb2;

end
