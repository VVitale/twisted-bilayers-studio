%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to choose more accurate Slater-Koster
% interlayer parameters for interlayer p-p and 
% pz-dz2 interactions. Interpolates to given 
% vertical distance for two atoms involved in
% hopping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% hop - Type of hopping parameter: 1: 'pp' or 2: 'pd'
% dsqr - Modulus of distance between 2 atoms
% rz - Magnitude of z-component of the connecting
%       vector between 2 atoms involved in hopping
% zpar - Interlayer seperations corresponding to each
%       parameter in intpar cell array
% intpar - SK parameter cell array for each 
%          interlayer seperation
%
% Outputs:
% Vsigval - Interpolated value of the SK parm. Vsig
% Vpival - Interpolated value of the SK param. Vpi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Kemal Atalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vsigval,Vpival] = interpolate_SK(hop,dsqr,rz,zpar,intpar)

% Range
zpar = cell2mat(zpar);

% Range of hopping values calculated at different parameters
Vsig_lst = zeros(size(intpar)); 
Vpi_lst = zeros(size(intpar)); 

% Calculate hoppings to put in hoprange
for ii = 1:size(intpar)
    % Choose the type of hopping values to interpolate
    if(hop == 1)
        % p-p SL parameters
        Vsig = intpar{ii}(1)*exp(-(dsqr/abs(intpar{ii}(3)))^intpar{ii}(5));
        Vpi = intpar{ii}(2)*exp(-(dsqr/abs(intpar{ii}(4)))^intpar{ii}(6));
    elseif(hop == 2)
        % pd SK parameters 
        Vsig = intpar{ii}(1)*(dsqr^intpar{ii}(3))*cos(intpar{ii}(5)*dsqr + intpar{ii}(7));
        Vpi  = intpar{ii}(2)*(dsqr^intpar{ii}(4))*cos(intpar{ii}(6)*dsqr + intpar{ii}(8));
    else
        error('Incorrect hopping type set in interpolate_SK.m')
    end
    Vsig_lst(ii) = Vsig; % add to array
    Vpi_lst(ii) = Vpi; % add to array
end

% Interpolate if it is in zpar range
% (might implement extrapolation in the future if out of range)
if(rz < zpar(1))
    Vsigval=Vsig_lst(1);
    Vpival=Vpi_lst(1);
elseif(rz > zpar(length(zpar)))
    Vsigval=Vsig_lst(length(zpar));
    Vpival=Vpi_lst(length(zpar));
else
    % Interpolate
    Vsigval=interp1(zpar,Vsig_lst,rz);
    Vpival=interp1(zpar,Vpi_lst,rz);
end

end
