 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to find a set of commensurate moire lattices 
% and the atomic coordinates for twisted homobilayers.
% The twisted bilayer is assumed to be made of 
% a material whose conventional cell is hexagonal, 
% e.g. graphene, hexagonal boron nitride (hBN) and 
% Transition-metal dichalcogenides (TMDs)
%
%
% INPUT: 
%  n      - Integer parameter that defines Moire cell,
%           m=n+1 by default, see Ref. 1
%  a      - lattice parameter of 2D material
%  c      - Interlayer distance when starting 2D material
%           consists of a monolayer, e.g. graphene, OR 
%           disance between central atoms of 2D material 
%           for trilayers, e.g. distance 
%           between transition metal atoms in the two
%           TMD layers.
%  comb   - Define which combination of layers to consider
%           Possible choices are: 'graphene', 'hBN', 
%           'MoS2', 'WSe2', 'MoSe2', 'WS2'
% plot_f  - Format for plotting ('xyz' or 'xsf')
% kpath   - Path in k-space. Implemented paths are:
%           * K-Gamma-M-K' -> 'KMGKp'
%           * M-Gamma-K-M  -> 'MGKM'
%           * Gamma-M-K-Gamma -> 'GMKG'
%           * Gamma-K-M    -> 'GKM'
%           * K-M          -> 'KM'
%           * Monkhorst-Pack -> 'uniform'
% init_nk - Total number of kpoints for uniform MP grids OR 
%           number of k-points in the first segment of the 
%           chosen path
%
% OUTPUT:
%  <rootname>.xyz(xsf)  - Cartesian coordinates in XYZ 
%                         (XSF) format
%  positions.<rootname>.dat  - Cartesian coordinates & 
%                              Moire lattice vectors
%                              to be used as input file
%                              for TB_free.x
%  lammps_positions.<rootname>.dat - Lammps geometry
%                                    file. 
%                                    N.B. The Moire 
%                                    cell is rotated
%                                    such as one of 
%                                    the vector is
%                                    parallel to
%                                    on the x-axis
%  potential.<rootname>.dat - Potential file to be
%                             used as input in 
%                             TB_free.x. Only for hBN
%                             at the moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale, October 2019
%
% Refs. 
%       1) Nano Lett. 10, 3, 804-808 (2010)
%       2) PRB 86, 155449 (2012) 
%
% version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 11-band TB model PRB 92 205108 (2015)
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      d(M-X)   d(X-X)    c (sqrt(d(M-X)^2-1/3*a^2)/(3/4-u))    u
% MoS2  3.1824  2.41      3.13    12.29
% WS2   3.1817  2.42      3.14    12.32
% MoSe2 3.3174  2.54      3.34    12.90
% WSe2  3.3155  2.55      3.35    12.96

% Set from SW classical potential LAMMPs
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      c       u 
% MoS2  3.110   12.08  0.621
% WS2   3.129   12.12  0.622
% MoSe2 3.3096  12.84  0.621
% WSe2  3.289   12.74  0.621
    
% Parameters from our DFT calculations
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          a      c         u
% MoS2  3.183  12.0264    0.621
% WS2   3.182  12.1634    0.621
% MoSe2 3.315  12.9567    0.621
% WSe2  3.3145 13.0089    0.621

%%%%%%%%%%%% BEGIN INPUT %%%%%%%%%%%%%
% Integer for constructing Moire supercell
n = 6;
% Lattice parameter of monolayer
a = 3.3145; %
% Interlayer distance
c = 13.0089/2;
% Chemical symbol for monolayer
comb = 'WSe2';
% Info for writing Cartesian coordinates and potential to file
write_positions = true;
write_potential = false;
write_cart      = false;
% Whether to XYZ or XSF format
plot_f = 'xyz';
write_lammps_input = false;
% Initial phase
phase = 'AA';

% Info for writing k-points and path in k-space
write_kpath = false;
kpath = 'GMKG';
init_nk = 10;
%%%%%%%%%%%% END INPUT %%%%%%%%%%%%%

switch comb
    case 'graphene'
        disp('Structure and input files for twisted bilayer graphene')
        disp(' ')
        [num,pos,pos2] = TBLG_cell(n,a,c,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'hBN'
        disp('Structure and input files for twisted bilayer hBN')
        disp(' ')
        [num,pos,pos2] = TBLhBN_cell(n,a,c,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'MoS2'
        % Atomic masses
        m1 = 95.94;
        m2 = 32.06;
        % Atoms info
        num_at1 = 42;
        num_at2 = 16;
        disp('Structure and input files for twisted bilayer MoS2')
        disp(' ')
        [num,pos,pos2] = TBLTMD_cell(n,a,c,'Mo','S',m1,m2,num_at1,num_at2,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'MoSe2'
        % Atomic masses
        m1 = 95.94;
        m2 = 78.96;
        % Atoms info
        num_at1 = 42;
        num_at2 = 34;
        disp('Structure and input files for twisted bilayer MoSe2')
        disp(' ')
        [num,pos,pos2] = TBLTMD_cell(n,a,c,'Mo','Se',m1,m2,num_at1,num_at2,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'WSe2'
        % Atomic masses
        m1 = 183.84;
        m2 = 78.96;
        num_at1 = 74;
        num_at2 = 34;
        disp('Generating structure and input files for twisted bilayer WSe2')
        disp(' ')
        [num,pos,pos2] = TBLTMD_cell(n,a,c,'W','Se',m1,m2,num_at1,num_at2,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'WS2'
        % Atomic masses
        m1 = 183.84;
        m2 = 32.06;
        num_at1 = 74;
        num_at2 = 16;
        disp('Generating structure and input files for twisted bilayer WSe2')
        disp(' ')
        [num,pos,pos2] = TBLTMD_cell(n,a,c,'W','S',m1,m2,num_at1,num_at2,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    otherwise
        error('Material not supported. Available choices are: graphene, hBN, MoS2 and WSe2')
end
%end
disp('Done! Have a nice day :-)')
