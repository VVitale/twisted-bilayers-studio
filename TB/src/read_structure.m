%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to read in geometry info from 
% structure.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% filename - name of geometry file ('structure.dat')
%
% Outputs:
% nat - Number of atoms on each layer
% theta - twist angle in degrees
% unit_cell - cell vectors of unit cell of monolayer (ang)
% cell - cell vectors of moire cell (ang)
% structure - structure containing info on all atoms
%             Atom label, layer, X, Y, Z (ang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nat,theta,strain,unit_cell,cell,structure,bondlength] = read_structure(filename,nlayer)
    fileID = fopen(filename,'r');
    data = textscan(fileID,'%f %f %f %f %f',1);
    % Read number of atoms
    for ilayer = 1 : nlayer
       nat(ilayer) = cell2mat(data(ilayer));
    end
    data = textscan(fileID,'%f %f %f %f %f',1);
    for ilayer = 1 : nlayer
       bondlength(ilayer) = cell2mat(data(ilayer));
    end
    % Read theta
    data = textscan(fileID,'%f %f %f %f %f',1);
    for ilayer = 1 : nlayer
       theta(ilayer) = cell2mat(data(ilayer));
    end
    % Read strain
    data = textscan(fileID,'%f %f %f %f',1);
    for ilayer = 1 : nlayer
       strain(ilayer) = cell2mat(data(ilayer));
    end
    bondlength = bondlength.*strain;
    % Read unit cell
    data = textscan(fileID,'%f %f %f',3);
    unit_cell = cell2mat(data);
    % Read cell
    data = textscan(fileID,'%f %f %f',3);
    cell = cell2mat(data);
    data = textscan(fileID,'%s %f %f %f %f',sum(nat));
    structure = cell2struct(data,{'name','layer','x','y','z'},2);
    fclose(fileID);    
end
