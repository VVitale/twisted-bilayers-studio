%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to find a set of commensurate moire lattices 
% and the atomic coordinates for two materials with 
% hexagonal cells and different lattice parameters.
%
%
% INPUT: 
%  a     - Lattice parameters of substrate
%  p     - Ratio between substrate's (bottom layer) 
%          lattice parameter and top layer's
%  theta - Target twist angle between layers in degrees
%  rt    - +/- range of twist angles around theta in degrees
%  rp    - Allowed error in p for the search
%  c      - Interlayer distance when starting 2D material
%           consists of a monolayer, e.g. graphene, OR 
%           disance between central atoms of 2D material 
%           for trilayers, e.g. distance 
%           between transition metal atoms in the two
%           TMD layers.
%  comb  - Define which combination of layers to consider
%          * 1 = Graphene on hBN
%          * 2 = Graphene on Pt(111)
%          * 3 = Graphene on MoS2 (not implemented yet)
%          * -1 = Twisted bilayer graphene
% plot_f - Format for plotting ('xyz' or 'xsf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale, October 2019
%
% Ref. 
%      1) New J. Phys. 16 (2014)
%
% version 1.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful parameters
% Conversion factor Degree to Rad
deg2rad = pi/180;
% Rad to Degree
rad2deg = 1.0/deg2rad;
% Maximum number of atoms in bilayer
nat_max = 100000;
% Maximum integer for m and n (where n > m+1)
max_int = 100;
% Step in twist angle
dt = 0.001;
% abs(x_m - x_n) < threshold
threshold = 5E-6;
% Tollerance for checking if atom belongs to moirè unit cell 
tol = 1*10^(-5);

%%%%%%%%%%%% BEGIN INPUT %%%%%%%%%%%%%
% lattice parameter of bottom layer in angstroms
a = 3.1824;
% lattice mismatch between bottom and top layer, i.e. a_bottom / a_top
p = 3.1824/3.1817;
% target twist angle in deg
theta = 4;
% range around target twist angle in deg
r_t = 1;
% range around target lattice ratio
r_p = 0.01;
% interlayer distance in angstroms (slightly different definition for TMDs, see HEADER)
c = 12.2985;
% used to define relative position of chalcogenides in TMDs
u = 0.621;
z = 0.0;
% number of atomic species
ntype = 4;
% Whether a layer contains two (graphene) or three (TMD) atoms per unit cell
layer1 = 'TMD';
layer2 = 'TMD';
phase = 'AA';

% Info for writing Cartesian coordinates and potential to file
write_positions = true;
write_potential = true;
write_cart      = true;
% Whether to XYZ or XSF format
plot_f = 'xsf';
write_lammps_input = true;

% Atomic symbol of each atomic species
id11 ="W";%'Mo';
id12 ="Se";%'S';
id13 ="Se";%'S';
id21 ="W";%'W';
id22 ="S";%'Se';
id23 ="S";%'Se';
% Atomic number of each species (used in XSF file)
at_num11 =74 ;%42;
at_num12 =34 ;%16;
at_num13 =34 ;%16;
at_num21 =74 ;%74;
at_num22 =16 ;%34;
at_num23 =16 ;%34;

% Mass of each atomic species in a.u. (used in LAMMPS input file)
mass_1 =183.84;%95.94;
mass_2 =78.96;%32.06;
mass_3 =183.84;%183.84;
mass_4 =32.06;%78.96;

% Potential for each atomic species
pot_11 = 0.0;
pot_12 = 0.0;
pot_13 = 0.0;
pot_21 = 0.0;
pot_22 = 0.0;
pot_23 = 0.0;

% For LRI only
%dmax = 4.50;
%dmin = 4.35;

% PATH in k-space
path = 'KGMKp';

%%%%%%%%%%%% END INPUT %%%%%%%%%%%%%


%%%%%%%%%%% START SEARCH %%%%%%%%%%


   if(layer1=="TMD")
      f1 = join([id11,id12,"2"],"");
   else
      f1 = join([id11,id12],"");
   end
   if(layer2=="TMD")
      f2 = join([id21,id22,"2"],"");
   else
      f2 = join([id21,id22],"");
   end
   
   if(layer1=="TMD" && layer2=="TMD")
      comb = 1;
   elseif(layer1=="graphene" && layer2=="graphene")
      comb = 2;
   elseif((layer1=="graphene" && layer2=="TMD") || (layer2=="graphene" && layer1=="TMD"))
      comb = 3;
   end

   filename = join([f1,f2],"_on_");
% Check inputs

% if (theta==0 && r_t==0)
%    m=64;
%    n=66;
%    [alatL, pos, pos2] = notwist(filename, 0.0, p, c, u, a, m, n, path);
% else
    size_t = 2*int64(r_t/dt) + 1;
    m_s = zeros(size_t,1);
    n_s = zeros(size_t,1);
    theta_s = zeros(size_t,1);
    p_s = zeros(size_t,1);
    if(phase=='AA')
       theta = theta;
    elseif(phase=='AB')
       theta = 60 - theta;
    end
    [m_s,n_s,p_s,theta_s,nmoire] = find_moire(theta,dt,r_t,max_int,size_t,1,threshold,p,r_p,m_s,n_s,p_s,theta_s);
% Check at least one structure has been found
if (nmoire == 0)
    msg = ['No commensurate moire superlattices could be found' ...
        ' for theta in range [',num2str(theta - r_t),';', ...
        num2str(theta + r_t),'].', ' Try with a different' ...
        ' initial twist angle or increase the range.'];
    disp(msg)
    error('Aborting')
end

nmoire_f = nmoire;
index = find(m_s);

% Print some useful info
disp(' ')
msg=['Total number of unique commensurate moirè lattice found: ', num2str(nmoire_f)];
disp(msg)
disp(' ')
% Generate unit cells and atomic coordinates for each moire superlattice
for it = 1:nmoire
    % Find index of non-zero elements of m_s    
    im = m_s(index(it));
    in = n_s(index(it));
    % Find number of atoms in moire cells
    corr_nat = num_atoms(int64(im),int64(in),int64(comb));
    msg = ['Generating Cartesian atomic coordinates for moirè lattice ',num2str(it),' ...'];
    disp(msg)
    % Check number of atoms is less than max allowed by TB 
    if (corr_nat < nat_max)
        itheta = theta_s(index(it));
        ip = p_s(index(it));
        st=sin(itheta);
        ct=cos(itheta);
        Rtheta = [ct -st; st ct];
        I = eye(2)*(1/ip);
        alat1 = [1 -0.5; 0 sqrt(3)/2]*a;
        alat2 = I*Rtheta*alat1;
        
        % Generate direct lattice
        intL = [im in;-in (im - in)];
        alatL = intL*alat1';
        ialatL = inv(alatL');
        
        % Useful for file names
        str_theta = num2str(itheta*rad2deg);
        str_m = num2str(im);
        str_n = num2str(in);

              
        % Check commensurability with top layer
        intL2 = alatL/(alat2');
        if (any(abs(intL2 - round(intL2))>0.001))
            error('Wrong moirè lattice vectors')
        end
        
        % Allocate and initialise arrays
        % +/- number of cells to loop on
        ncell=100; 
        % Total number of neighboring cells
        neigh=2*ncell + 1; 
        % xyz coordinates of bottom layer
        pos=zeros(corr_nat,3); 
        % xyz coordinates of top layer
        pos2=zeros(corr_nat,3); 
        pos(:,:) = -99999; 
        pos2(:,:) = -99999;
        % Array containing the atoms' labels for bottom and top layers
        id1(1,1:corr_nat) = "X"; 
        id2(1,1:corr_nat) = "X";
        % Array containing the atoms' numbers
        at_num1(1,1:corr_nat) = 0;
        at_num2(1,1:corr_nat) = 0;

        % CYCLE OVER NEIGHBORING CELLS
        % Bottom layer first
        [nat,nat1,nat2,pos,pos2,id1,id2,at_num1,at_num2]=generate_structure(a,u,c,z,I,Rtheta,ncell,alat1,ialatL,tol,pos,pos2,id1,id2,id11,id12,id13,id21,id22,id23,at_num1,at_num2,at_num11,at_num12,at_num13,at_num21,at_num22,at_num23,layer1,layer2,phase);

        disp('done')
        % Check number of atoms is correct
        if (nat ~= corr_nat)
            msg = ['Number of atoms is incorrect. Expected ',num2str(corr_nat),' found ',num2str(nat)];
            disp(msg)
            error('Error')
        else
            msg = ['Number of atoms is correct : ', num2str(nat)];
            disp(msg)
        end
       

        if (write_kpath)
           msg = ['Generating kpoint grid for moirè lattice ',num2str(it),' ...'];
           disp(msg)
           % Generate reciprocal lattice
           recL = 2*pi*inv(alatL');
           irecL = inv(recL);
           trecL = recL'; 
           ralat1 = 2*pi*inv(alat1');
           ralat2 = 2*pi*inv(alat2');
           write_kpts(15,ralat2,recL,str_theta,str_m,str_n,path,2)
           disp('done')
        end
        
        
        pos = pos(1:nat1,:);
        pos2 = pos2(1:nat2,:);
        id1 = id1(1,1:nat1);
        id2 = id2(1,1:nat2);

        % LRI
        %rlri1 = zeros(nat1,1);
        %rlri2 = zeros(nat2,1);
        %[rlri1,rlri2] = compLRIghBN(pos,pos2,nat1,nat2,alatL,a,id1,ip,false,1);
        %delta = 0.5*(dmax-dmin);
        %d1 = zeros(nat1,1);
        %d2 = zeros(nat2,1);


        %%%% WRITE POSITION FILE FOR TB_FREE.X %%%%
        if (write_positions)
           write_pos(nat1,nat2,filename,alatL,theta*180/pi,alat1,id1,id2,pos,pos2,str_m,str_n,str_theta)
        end

        %%%% WRITE POTENTIAL FILE FOR TB_FREE.X %%%%
        if (write_potential)
           write_pot(nat1,nat2,filename,id1,id2,id11,id12,id13,id21,id22,id23,pot_11,pot_12,pot_13,pot_21,pot_22,pot_23,pos,pos2,str_m,str_n,str_theta)
        end
             
        %%%% WRITE COORDINATES FILE %%%%
        if(write_cart)
           if (plot_f=='xyz')
              write_xyz(nat1,nat2,filename,id1,id2,pos,pos2,str_m,str_n,str_theta)
           elseif (plot_f=='xsf')
              write_xsf(nat1,nat2,filename,alatL,c,at_num1,at_num2,pos,pos2,str_m,str_n,str_theta)
           end
        end

        %%%% WRITE LAMMPS INPUT FILE %%%%
        if (write_lammps_input)
           m(1) = mass_1;
           m(2) = mass_2;
           m(3) = mass_3;
           m(4) = mass_4;

           write_lammps(nat1,nat2,filename,alatL,c,ntype,m,id1,id2,id11,id12,id13,id21,id22,id23,pos,pos2,str_m,str_n,str_theta)
        end

        %figure
    else
        nmoire_f = nmoire_f - 1;
        index(it) = 0;
    end
end


function y = num_atoms(m,n,c)
    if ~isinteger(m)
        error('First imput must be an integer (m)')
    end
    if ~isinteger(n)
        error('Second imput must be an integer (n)')
    end
    if ~isinteger(c)
        error('Third imput must be an integer (c)')
    elseif c< -1 || c > 3 || c==0
        error('c must be equal to: 1 or 2 or 3 or -1')
    end
    
    if (c==1) 
        y = 6*(n^2 - n*m + m^2) + 3*(1 - n + 2*m);
    elseif (c==2)
        y = 4*(n^2 - n*m + m^2) + 2*(1 - n + 2*m);
    elseif (c==3)
        y = 3*(n^2 - n*m + m^2) + 3/2*(1 - n + 2*m) + 2*(n^2 -n*m + m^2) + (1 -n + 2*m) - 3;
    end
    
end
