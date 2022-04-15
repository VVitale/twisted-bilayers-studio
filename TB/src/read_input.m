%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to read in input file input.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% filename - name of input file ('input.dat')
%
% Outputs:
% - task                   (I) 1: bandstructure  2: Chern number and Berry curvature   3: DOS  4: charge density   5: exciton spectra and JDOS
%                              6: Plot wavefunctions at Gamma, M and K   7: Build Hamiltonian on uniform grid and write in hdf5 format
% - restart                (L) Whether this is a restart calculation
% - read_ham               (L) Whether to read Hamiltonian 
% - geomfname              (S) Name of geometry file
% - gw_par                 (L) Whether to use GW parameters
% - d12                    (R) maximum distance between chalcogen atom on bottom layer and chalcogen atom on top layer for computing 
%                              interlayer hopping (aXX2 = 5.0)
% - nlayer                 (I) Number of layers
% - tmdc                   (I) TMD material code --- 1: MoS2   2: MoSe2   3: WS2   4:WSe2, e.g. MoS2 monolayer tmdc: = 1, MoS2 bilayer: tmdc = 11, 
%                              MoS2 trilayer tmdc : 111, WS2/MoS2 bilayer: tmdc = 13
% - convention             (I) Convention for the unit cell vectors of the monolayer (a=sqrt(3)*aXM)
%                                                          a1                   a2
%                              convention 1 (120 deg) | (1, 0)*a   |  (-1/2,sqrt(3)/2)*a
%                              convention 2 (60 deg)  | (1, 0)*a   |  (1/2,sqrt(3)/2)*a
% - knum                   (I) Number of k-points in the first segment (Gamma-M) when task = 1. Number of k-points along one reciprocal
%                              vector when task = {2,3,4,5,7}, in this case the total number of k-points is knum^2
% - miniBZ                 (L) Whether to use Gamma-M-K-Gamma of the miniBZ of Gamma-M-K-Gamma of the BZ of the monolayer folded onto the miniBZ
% - interlayer_int         (L) Whether to compute the interlayer interaction
% - spin_orbit             (L) Whether to compute spin-orbit 
% - reduced_workspace      (L) Whether to use eig or eigs in the diagonalisation
% - num_eigs               (I) Number of conduction bands when reduced_workspace = true 
% - eigvecs                (L) Whether to compute the eigenvectors
% - save_workspace         (L) Whether to save the orbital class for restart calculations
% - num_workers            (I,S) Number of workers in parallel diagonalisation or 'max' maximum number of workers allowed by local profile
% - dos_spacing            (R) Energy spacing in eV for DOS plot (e.g. 0.005)
% - gradient               (L) Whether to compute the x component of the Hamiltonian gradient
% - Interpd                (L) Whether to compute the pz-dz2 hopping
% - lambda                 (R) Shift of eigenspectrum, H = H + lambda*I
% - bse_num_vbnd           (I) Number of valence bands in BSE
% - bse_num_cbnd           (I) Number of conduction bands in BSE
% - bse_int                (L) Whether to compute spectra with BSE interaction
% - bse_broadening         (R) Broadening of the Lorentzian used to approximate the delta function in optical spectra
% - bse_eig_plot           (I) Number of excitonic eigenvfunction to plot
% - bse_irh                (I) Index of hole atomic position in the (moire) unit cell
% - bse_serial             (L) Whether to run the BSE solver serially
% - bse_shifted            (L) Whether to use a randomly shifted from (0,0) uniform k-grid for BSE
% - write_ham              (L) Write Hamiltonian to hdf5 format
% - ham_fname              (S) Name of hamiltonian hdf5 file
% - flipped                (S) Whether a layer has been flipped to form a 2H stacking
% - read_kpts              (L) Whether to read k-points
% - ef_strength            (R) Electric-field strength in V/Ang

% R = real number
% I = integer number
% L = boolean (true, false)
% S = string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, restart, read_ham, geomfname, GW, aXX2, ...
    nlayer, tmdc, convention, knum, miniBZ, interlayer_int, Interpd, lsoc, ...
    reduced_workspace, num_cond, compute_eigvecs, save_workspace, ...
    num_workers, dE, gradient, lambda, bse_num_vbnd, bse_num_cbnd, bse_int, ...
    bse_broadening, bse_eig_plot, bse_irh, bse_serial, bse_shifted, ...
    flipped, write_ham, ham_fname, read_kpts, ef_strength] = read_input(filename)

    keys = {'task';                  %1
            'geom_file';             %2
            'd12';                   %3
            'nlayer';                %4
            'tmdc';                  %5
            'convention';            %6
            'knum';                  %7
            'interlayer_int';        %8
            'num_workers';           %9
            'restart';               %10
            'read_ham';              %11
            'gw_par';                %12
            'mini_bz';               %13
            'spin_orbit';            %14
            'reduced_workspace';     %15
            'num_eigs';              %16
            'eigvecs';               %17
            'save_workspace';        %18
            'dos_spacing';           %19
            'gradient';              %20
            'Interpd';               %21
            'lambda';                %22
            'bse_num_vbnd';          %23
            'bse_num_cbnd';          %24
            'bse_int';               %25
            'bse_broadening';        %26
            'write_ham';             %27
            'ham_fname';             %28
            'flipped';               %29
            'read_kpts';             %30
            'bse_eig_plot';          %31
            'bse_irh';               %32
            'bse_serial';            %33
            'bse_shifted';           %34
            'ef_strength';};         %35
    numeric_keys = [1,3,4,5,6,7,9,16,19,22,23,24,26,29,31,32,35];
    fundamental_keys = [1,2,4,5,6,7,9,29];
    len_keys = size(keys,1);
    value = cell(len_keys,1);

    % Set default values
    default_values = {1;                         %1
                      'positions_example.dat';   %2
                      5;                         %3
                      1;                         %4
                      [];                        %5
                      1;                         %6
                      0;                         %7
                      'false';                   %8
                      1;                         %9
                      'false';                   %10
                      'false';                   %11
                      'false';                   %12
                      'false';                   %13
                      'false';                   %14
                      'false';                   %15
                      0;                         %16
                      'false';                   %17
                      'false';                   %18
                      0.05;                      %19
                      'false';                   %20
                      'false';                   %21
                      0;                         %22
                      0;                         %23
                      0;                         %24
                      'false';                   %25
                      0.1;                       %26
                      'false';                   %27
                      'ham.h5';                  %28
                      [];                        %29
                      'false';                   %30
                      1;                         %31
                      1;                         %32
                      'false';                   %33
                      'false';                   %34
                      0.0;};                     %35

    loc_found = 0;
    i = 0;
    found = false;
    %fileID = fopen(filename,'r');
    input = textread(filename,'%s','delimiter','\n');
    nlines = size(input,1);
    for iline = 1 : nlines
        if(~isempty(input{iline}))
           tmp = textscan(input{iline},'%s%f%s', 'Delimiter',':!', 'MultipleDelimsAsOne',true);
           tmp_key = deblank(tmp{1}{1});
           % Check spelling is correct
           [found,id] = ismember(tmp_key,keys);
           if(found)
              i = i + 1;
              loc_found(i) = id;
              if(any(numeric_keys == id))
                 tmp = textscan(input{iline},'%s%f%s', 'Delimiter',':!', 'MultipleDelimsAsOne',true);
              else
                 tmp = textscan(input{iline},'%s%s%s', 'Delimiter',':!', 'MultipleDelimsAsOne',true);
              end
              value{id} = tmp{2};
           else
              disp(['Error in read_input.m: wrong keyword ',tmp_key])
              error('Aborting ...')
           end
        end
    end

    nnz_mem = ismember(fundamental_keys,loc_found);
    check_nnz = nnz(nnz_mem);
    if(check_nnz < length(fundamental_keys))
        for ii = 1 : length(nnz_mem)
           fprintf('Cannot find "%s" keyword !', keys{ii})
        end
        error('Aborting ...')
    end

   for ikey = 1 : length(keys)
       if(~ismember(ikey,loc_found))
           value{ikey} = default_values{ikey};
       end
   end

    % Set input values
    task              = value{1};
    geomfname         = deblank(char(value{2}));
    aXX2              = value{3};
    nlayer            = value{4};
    tmp               = num2str(value{5});
    tmdc = zeros(nlayer,1);
    for il = 1 : nlayer
    	tmdc(il)    = str2num(tmp(il));
    end
    convention        = value{6};
    knum              = value{7};
    interlayer_int    = strcmp(deblank(char(value{8})),'true');
    num_workers       = value{9};
    restart           = strcmp(deblank(char(value{10})),'true');
    read_ham          = strcmp(deblank(char(value{11})),'true');
    GW                = strcmp(deblank(char(value{12})),'true');
    miniBZ            = strcmp(deblank(char(value{13})),'true');
    lsoc              = strcmp(deblank(char(value{14})),'true');
    reduced_workspace = strcmp(deblank(char(value{15})),'true');
    num_cond          = value{16};
    compute_eigvecs   = strcmp(deblank(char(value{17})),'true');
    save_workspace    = strcmp(deblank(char(value{18})),'true');
    dE                = value{19};
    gradient          = strcmp(deblank(char(value{20})),'true');
    Interpd           = strcmp(deblank(char(value{21})),'true');
    lambda            = value{22};
    bse_num_vbnd      = value{23};
    bse_num_cbnd      = value{24};
    bse_int           = strcmp(deblank(char(value{25})),'true');
    bse_broadening    = value{26};
    write_ham         = strcmp(deblank(char(value{27})),'true');
    ham_fname         = deblank(char(value{28}));
    clear tmp
    format            = join(['%0',num2str(nlayer),'i']);
    tmp               = num2str(value{29},format);
    flipped = zeros(nlayer,1);
    for il = 1 : nlayer
    	flipped(il)    = str2num(tmp(il));
    end
    read_kpts         = strcmp(deblank(char(value{30})),'true');
    bse_eig_plot      = value{31};
    bse_irh           = value{32};
    bse_serial        = strcmp(deblank(char(value{33})),'true');
    bse_shifted       = strcmp(deblank(char(value{34})),'true');
    ef_strength       = value{35};

    % Check input
    if(~isnumeric(aXX2))
       error('Error in read_input.m: d12 must be a real number. Aborting ...')
    else
       if(aXX2<=0)
          error('Error in read_input.m: d12 must be a real number greater than 0.0. Aborting ...')
       end
    end
    
    if(~isnumeric(dE))
       error('Error in read_input.m: dos_spacing must be a real number. Aborting ...')
    else
       if(dE<=0)
          error('Error in read_input.m: dos_spacing must be a real number greater than 0.0. Aborting ...')
       end
    end
    
    if(~isnumeric(lambda))
       error('Error in read_input.m: lambda must be a real number. Aborting ...')
    end

    if(~isnumeric(bse_broadening))
       error('Error in read_input.m: bse_broadening must be a real number. Aborting ...')
    end

    % INTEGER INPUTS
    if(mod(convention,1)~=0)
       error('Error in read_input.m: convention must be an integer. Aborting ...')
    else
       if(convention~=1 && convention~=2)
          error('Error in read_input.m: convention must be either 1 or 2. Aborting ...')
       end
    end
    
    if(mod(nlayer,1)~=0)
       error('Error in read_input.m: tmd1 must be an integer. Aborting ...')
    else
       %if(tmdc1~=1 && tmdc1~=2 && tmdc1~=3 && tmdc1~=4)
       %   error('Error in read_input.m: tmd1 must be an integer between 1 and 4. Aborting ...')
       %end
    end
    
    %if(mod(tmdc,1)~=0)
    %   error('Error in read_input.m: tmd2 must be an integer. Aborting ...')
    %else
    %   if(tmdc2~=1 && tmdc2~=2 && tmdc2~=3 && tmdc2~=4)
    %      error('Error in read_input.m: tmd2 must be an integer between 1 and 4. Aborting ...')
    %   end
    %end
    
    if(mod(task,1)~=0)
       error('Error in read_input.m: task must be an integer. Aborting ...')
    else
       if(task~=1 && task~=2 && task~=3 && task~=4 && task~=5 && task~=6 && task~=7)
          error('Error in read_input.m: task must be an integer between 1 and 7. Aborting ...')
       end
    end
    
    if(mod(knum,1)~=0)
       error('Error in read_input.m: convention must be an integer. Aborting ...')
    end
    
    if(mod(num_cond,1)~=0)
       error('Error in read_input.m: num_eigs must be an integer. Aborting ...')
    end
    
    if(mod(num_workers,1)~=0)
       error('Error in read_input.m: num_workers must be an integer or max. Aborting ...')
    else
       if(num_workers==-1)
          num_workers = 'max';
       end
    end
    
    if(mod(bse_num_vbnd,1)~=0)
       error('Error in read_input.m: bse_num_vbnd must be an integer. Aborting ...')
    end  

    if(mod(bse_num_cbnd,1)~=0)
       error('Error in read_input.m: bse_num_cbnd must be an integer. Aborting ...')
    end
  
    % LOGICAL INPUTS
    if(~islogical(interlayer_int))
       error('Error in read_input.m: interlayer_int must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(lsoc))
       error('Error in read_input.m: spin-orbit must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(miniBZ))
       error('Error in read_input.m: mini_bz must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(reduced_workspace))
       error('Error in read_input.m: reduced_workspace must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(compute_eigvecs))
       error('Error in read_input.m: eigvecs must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(save_workspace))
       error('Error in read_input.m: save_workspace must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(restart))
       error('Error in read_input.m: restart must be a logical variable, either false or true. Aborting ...')
    end
    
    if(~islogical(read_ham))
       error('Error in read_input.m: read_ham must be a logical variable, either false or true. Aborting ...')
    end

    if (~islogical(Interpd))
       error('Error in read_input.m: Interpd must be a logical variable, either false or true. Aborting ...')
    end

    if (~islogical(bse_int))
       error('Error in read_input.m: bse_int must be a logical variable, either false or true. Aborting ...')
    end

    for il = 1 : nlayer
       if (flipped(il) ~=0 && flipped(il) ~=1)
          error('Error in read_input.m: flipped must be either 0 or 1. Aborting ...')
       end
    end

    % Check input
    % Check if geometry file exists
    if (~exist(geomfname,'file'))
       error('Cannot find geometry file %s. Aborting ...', geomfname)
    end

    if(~restart && read_ham)
       error('Error in read_input.m: read_ham cannot be true if restart is false. Aborting ...')    
    end 

    % Check compute_eigvecs=true when computing charge density
    if((task==2 || task==3 || task==4 || task==5 || task==6) && ~compute_eigvecs)
      fprintf('WARNING: task %i required but compute_eigvecs = false. Setting compute_eigvecs = true.',task)
      compute_eigvecs = true;
    end

    if(task == 5)
       gradient = true;
       if(mod(knum,2)==0)
          fprintf('WARNING: task %i requires knum to be and odd integer. Setting knum=knum + 1.', task)
          knum = knum + 1;
       end
    end

    if(task == 7 && ~write_ham)
      fprintf('WARNING: task %i required but write_ham = false.',task)
    end

    clear input tmp loc_found fundamental_keys numeric_keys keys default_values value;
end
