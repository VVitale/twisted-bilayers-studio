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
% WSe2 from PRB 60, 20 (1999) Kruger et al.
% a = 3.285
% c = 12.748
% MoS2 from PRB 84, 155413 (2011) Molina-Sanchez et al.
% a = 3.126
% c = 12.066
% WS2
% a = 3.125
% c = 12.145
% MoSe2 from RSC Adv. 6 5767 (2016) Peng et al.
% a = 3.300
% c = -

% Set from Ann. Phys. 526, 9-10, 347-357 (2014) Roldan et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      u'        c'    c (2c')    u (3/4 - u'/c)
% MoS2  3.160  1.586   6.140    12.280     0.621
% WS2   3.153  1.571   6.160    12.320     0.622
% MoSe2 3.288  1.664   6.451    12.902     0.621
% WSe2  3.260  1.657   6.422    12.844     0.621

% Experimental from PRB 89, 075409 (2014) He et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          a      c
% MoS2  3.160  12.29
% WS2   3.153  12.32
% MoSe2 3.288  12.90
% WSe2  3.282  12.96


% Set from  Phys. Chem. Chem. Phys 19, 13333 (2017) Morales Garcia et al. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      d(M-X)   d(X-X)    c (sqrt(d(M-X)^2-1/3*a^2)/(3/4-u))    u 
% MoS2  3.183   2.405    3.151   12.0264                               0.621
% WS2   3.182   2.416    3.139   12.1634                               0.622
% MoSe2 3.315   2.541    3.341   12.9567                               0.621
% WSe2  3.314   2.545    3.355   13.0089                               0.621

% Set from SW classical potential LAMMPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      c (sqrt(d(M-X)^2-1/3*a^2)/(3/4-u))    u 
% MoS2  3.110   12.08                                 0.621
% WS2   3.129   12.12                                 0.622
% MoSe2 3.3096  12.84                                 0.621
% WSe2  3.289   12.74                                 0.621

% Parameters from our DFT calculations
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          a      c         u
% MoS2  3.183  12.0264    0.621
% WS2   3.182  12.1634    0.621
% MoSe2 3.315  12.9567    0.621
% WSe2  3.3145 13.0089    0.621

% lattice parameter of bottom layer in angstroms
a = 3.183;
% lattice mismatch between bottom and top layer, i.e. a_bottom / a_top
p = 3.183/3.315;
% target twist angle in deg
theta = 3.389;
% range around target twist angle in deg
r_t = 0;
% range around target lattice ratio
r_p = 0.04*a;
% interlayer distance in angstroms (slightly different definition for TMDs, see HEADER)
c1 = 12.0264;
c2 = 12.9567;
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
write_potential = false;
write_cart      = true;
% Whether to XYZ or XSF format
plot_f = 'xsf';
write_lammps_input = true;
write_onetep_input = false;
mode = 'bs';
% PATH in k-space
write_kpath = true;
path = 'GMKG';

% Atomic symbol of each atomic species
id11 ="Mo";%'Mo';
id12 ="S";%'S';
id13 ="S";%'S';
id21 ="Mo";%'W';
id22 ="Se";%'Se';
id23 ="Se";%'Se';

% Atomic number of each species (used in XSF file)
if (id11 == 'Mo')
    at_num11 =42 ;%42;
    m(1) = 95.94;
elseif(id11 == 'W')
    at_num11 = 74;
    m(1) = 183.84;

end
if (id12 == 'S')
    at_num12 = 16;
    at_num13 = 16;
    m(2) = 32.06;
elseif (id12 == 'Se')
    at_num12 = 34;
    at_num13 = 34;
    m(2) = 78.96;
end
if(id21 == 'Mo')
    at_num21 = 42;
    m(3) = 95.94;
elseif(id21 == 'W')
    at_num21 = 74;
    m(3) = 183.84;
end

if(id22 == 'S')
    at_num22 =16 ;
    at_num23 =16 ;
    m(4) = 32.06;
elseif(id22 == 'Se')
    at_num22 = 34 ;%34;
    at_num23 = 34 ;%34;
    m(4) = 78.96;
end

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

   filename = join([f2,f1],"_on_");
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
        [nat,nat1,nat2,pos,pos2,id1,id2,at_num1,at_num2]=generate_structure(a,u,c1,c2,z,I,Rtheta,ncell,alat1,ialatL,tol,pos,pos2,id1,id2,id11,id12,id13,id21,id22,id23,at_num1,at_num2,at_num11,at_num12,at_num13,at_num21,at_num22,at_num23,layer1,layer2,phase);

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
           %recL = 2*pi*inv(alatL')
           %irecL = inv(recL);
           %trecL = recL'; 
           %ralat1 = 2*pi*inv(alat1');
           %ralat2 = 2*pi*inv(alat2');
           %write_kpts(15,ralat2,recL,str_theta,str_m,str_n,path,2)
           a1 = [alat1(:,1)' 0];
           a2 = [alat1(:,2)' 0];
           a3 = [0 0 40];
           ma1 = [alatL(1,:) 0];
           ma2 = [alatL(2,:) 0];
           ma3 = [0 0 40];
           v=abs(dot(a1,cross(a2,a3)));
           b1=2*pi*cross(a2,a3)/v;
           b2=2*pi*cross(a3,a1)/v;
           b3=2*pi*cross(a1,a2)/v;
           % Reciprocal lattice vectors of moire cell
           mv=abs(dot(ma1,cross(ma2,ma3)));
           mb1=2*pi*cross(ma2,ma3)/mv;
           mb2=2*pi*cross(ma3,ma1)/mv;
           mb3=2*pi*cross(ma1,ma2)/mv;
           recL=[mb1;mb2;mb3]
           k1=[0,0,0];
           k2=b1/2;
           k3=(2*b1-b2)/3;
           fM = k2*inv(recL);
           fK = k3*inv(recL);
           k2 = fM*recL
           C3 = [cos(pi/3) sin(pi/3) 0; -sin(pi/3) cos(pi/3) 0; 0 0 1];
           %C32 = [cos(2*pi/3) sin(2*pi/3) 0; -sin(2*pi/3) cos(2*pi/3) 0; 0 0 1];
           k3_all = [fK*recL;fK*recL*C3;fK*recL*C3*C3;fK*recL*C3*C3*C3;fK*recL*C3*C3*C3*C3;fK*recL*C3*C3*C3*C3*C3]
           [minv,iv] = min(sqrt(sum((k3_all - k2).^2,2)))
           k3 = k3_all(iv,:)
           fK = k3*inv(recL)
           if(abs(fM(1)) >= 1)
               if(fM(1)  > 0)
                   fM(1) = fM(1) - floor(fM(1)); %
               else
                   fM(1) = fM(1) - ceil(fM(1)); %
               end
           end
           if(abs(fM(2)) >= 1)
               if(fM(2)  > 0)
                   fM(2) = fM(2) - floor(fM(2)); %
               else
                   fM(2) = fM(2) - ceil(fM(2)); %
               end
           end
           if(abs(fK(1)) >= 1)
               if(fK(1)  > 0)
                   fK(1) = fK(1) - floor(fK(1)); %
               else
                   fK(1) = fK(1) - ceil(fK(1)); %
               end
           end
           if(abs(fK(2)) >= 1)
               if(fK(2)  > 0)
                   fK(2) = fK(2) - floor(fK(2)); %
               else
                   fK(2) = fK(2) - ceil(fK(2)); %
               end
           end
           fM
           fK
           k2 = fM*recL;
           k3 = fK*recL;
           disp('done')
        end
        
        
        pos = pos(1:nat1,:);
        pos2 = pos2(1:nat2,:);
        id1 = id1(1,1:nat1);
        id2 = id2(1,1:nat2);
        strain = p/ip;

        % LRI
        %rlri1 = zeros(nat1,1);
        %rlri2 = zeros(nat2,1);
        %[rlri1,rlri2] = compLRIghBN(pos,pos2,nat1,nat2,alatL,a,id1,ip,false,1);
        %delta = 0.5*(dmax-dmin);
        %d1 = zeros(nat1,1);
        %d2 = zeros(nat2,1);

        if(phase=='AB')
           itheta = pi/3 - itheta;
           str_theta = num2str(itheta*rad2deg);
       end

        %%%% WRITE POSITION FILE FOR TB_FREE.X %%%%
        if (write_positions)
           write_pos(nat1,nat2,filename,alatL,itheta*rad2deg,strain,alat1',id1,id2,pos,pos2,str_m,str_n,str_theta)
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
              %write_xsf(nat1,nat2,filename,rot_alatL,c1,at_num1,at_num2,pos,pos2,str_m,str_n,str_theta)
              write_xsf(nat1,nat2,filename,alatL,c1,at_num1,at_num2,pos,pos2,str_m,str_n,str_theta)
           end
        end

        phi=acos(alatL(1,1)/norm(alatL(1,:)));
        cf=cos(-phi);
        sf=sin(-phi);
        Rphi=[cf -sf; sf cf];
        rot_alatL = alatL*Rphi';
        rot_alat1 = alat1*Rphi';
        rot_pos = pos(:,1:2)*Rphi';
        rot_pos2 = pos2(:,1:2)*Rphi';
        pos = [rot_pos(1:nat1,:) pos(:,3)];
        pos2 = [rot_pos2(1:nat2,:) pos2(:,3)];

        %%%% WRITE LAMMPS INPUT FILE %%%%
        if (write_lammps_input)
           write_lammps(nat1,nat2,filename,rot_alatL,c1,ntype,m,id1,id2,id11,id12,id13,id21,id22,id23,pos,pos2,str_m,str_n,str_theta,phi)
        end

        if (write_onetep_input)
           write_onetep(nat1,nat2,filename,rot_alatL,m,id1,id2,id11,id12,id21,id22,pos,pos2,str_m,str_n,str_theta,mode)
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
