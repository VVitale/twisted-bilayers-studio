%clear all
%clc
fprintf('|-----------------------------------------------|\n')
fprintf('| Tight-binding code for twisted multilayer     |\n')
fprintf('| transition metal dichalcogenides, based on    |\n')
fprintf('| the original paper by Fang et al. [1]         |\n')
fprintf('|                                               |\n')
fprintf('| The code computes the bandstructure and the   |\n')
fprintf('| the DOS of both monolayer and twisted and     |\n')
fprintf('| untwisted multilayers.                        |\n')
fprintf('| Parameters obtained by Kemal Atlar using      |\n')
fprintf('| Quantum ESPRESSO and Wannier90, following     |\n')
fprintf('| the recipe outlined in Ref. 1                 |\n')
fprintf('|-----------------------------------------------|\n')
fprintf('| Inputs:                                       |\n')
fprintf('|  - input.dat Input file                       |\n')
fprintf('|  - geometry file containing the the structural|\n') 
fprintf('|    data of the system, name is defined in     |\n')
fprintf('|    input.dat file                             |\n')
fprintf('|-----------------------------------------------|\n')
fprintf('| Timeline                                      |\n')
fprintf('| February 2020 v1.0                            |\n')
fprintf('| March 2020 Added heterobilayer                |\n')
fprintf('| March 2020 Added relaxed geometries           |\n')
fprintf('| March 2020 Added SOC                          |\n')
fprintf('| March 2020 Added plot of the density          |\n')
fprintf('| April 2020 v1.1                               |\n')
fprintf('| April 2020 Added Berry curvature and Chern no.|\n')
fprintf('| April 2020 Excitonic spectra (BSE) and JDOS   |\n')
fprintf('| November 2022 v2.0                            |\n')
fprintf('|-----------------------------------------------|\n')
fprintf('| [1] Fang et al. PRB 92, 205108 (2015)         |\n')
fprintf('|-----------------------------------------------|\n')
fprintf('|                                               |\n')
fprintf('| Written by Valerio Vitale and Kemal Atalar    |\n')
fprintf('|                                               |\n')
fprintf('|-----------------------------------------------|\n\n')

t = clock;
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
fprintf('%i %s %i %02i:%02i\n\n',t(3),months{t(2)},t(1),t(4),t(5))
tic

%%%%%%%% INPUT %%%%%%%%%%%%%
restart = false;
if(exist(inpfname,'file'))
   fprintf('--> Reading input from input.dat file ... ')
   [task, restart, read_ham, geomfname, gw_par, aXX2, ...
    nlayer,tmdc, convention, knum, miniBZ, interlayer_int, Interpd, ...
    lsoc, reduced_workspace, num_cond, compute_eigvecs, save_workspace, ...
    num_workers, dE, gradient, lambda, vn, cn, elho_int, dL, bse_eig_plot, ...
    bse_irh, bse_serial, bse_shifted, flipped, ...
    write_ham, ham_fname, read_kpts, ef_strength]  = read_input(inpfname);
   fprintf('done\n\n')
else
   error('Input file not found! Aborting ...')
end
%%%%%%%% END %%%%%%%%%%%%%%


if (~onlydir)
fprintf('--> Checking number of workers ... ')
% Suppress output from parpool
pc = parcluster('local'); % If no pool, do not create new one.
if isempty(pc)
   poolsize = 0;
 else
   poolsize = pc.NumWorkers;
end
 
% Check number of workers is less than allowed number of workers set in 'local' profile
if (isnumeric(num_workers))
   if (num_workers > poolsize)
      fprintf('num_workers grater than allowed by local profile (%i)', poolsize)
      error('Aborting ...')
   end
else
   if(strcmp(num_workers,'max'))
      num_workers = poolsize;
   end
end

fprintf('done. Running on %i workers\n\n',num_workers)
end

% Check if running on a cluster with PBS queueing system
on_cluster = false;
switch getenv('PBS_ENVIRONMENT')

case 'PBS_BATCH' % PBS batch
   % Save the results
   filepath = strcat(getenv('TMPDIR'));%,getenv('PBS_JOBID'),'.save');
   % explicitly set the JobStorageLocation to the temp directory that was
   % created in your sbatch script
   pc.JobStorageLocation = strcat(filepath)
   on_cluster = true
end

parallel = false;
if(num_workers > 1 || on_cluster)
   parallel = true;
end

outfname = '';
for ilayer = 1 : nlayer
   switch tmdc(ilayer)
       case 1
           TMD{ilayer} = 'MoS2';
           metal{ilayer} = 'Mo';
           chalc{ilayer} = 'S';
       case 2
           TMD{ilayer} = 'MoSe2';
           metal{ilayer} = 'Mo';
           chalc{ilayer} = 'Se';
       case 3
           TMD{ilayer} = 'WS2';
           metal{ilayer} = 'W';
           chalc{ilayer} = 'S';
       case 4
           TMD{ilayer} = 'WSe2';
           metal{ilayer} = 'W';
           chalc{ilayer} = 'Se';
   end
outfname = join([outfname,TMD{ilayer},'_']);
end


if(~restart)
   fprintf('*** Starting calculation ***\n\n')
else
   fprintf('*** Starting a restart calculation ***\n\n')
end

fprintf('--> Reading geometry from %s file ... ', geomfname)
% Read in structure, number of atoms and lattice parameters from geomfname
[natoms,theta,strain,unit_cell,mcell,structure,bondlength] = read_structure(geomfname,nlayer);
fprintf('done\n\n')
tot_natoms = sum(natoms);
bondlength = bondlength./sqrt(3);
rel_pos = [structure.x,structure.y,structure.z];
% Read in structure, number of atoms and lattice parameters from pristine geomfname
geomfname_pristine = join([extractBefore(geomfname,'.dat'),'_pristine.dat']);

if(exist(geomfname_pristine,'file'))
   [rel_pos_ws,rel_id_ws,IN_ws,tr_vec] = wigner_seitz_cell(convention,structure,mcell,tot_natoms);
   rel_pos = rel_pos - tr_vec;
   [natoms_pristine,~,~,~,~,structure_pristine, ~] = read_structure(geomfname_pristine,nlayer);
   pris_pos = [structure_pristine.x,structure_pristine.y,structure_pristine.z];
   for iat = 1 : tot_natoms
      if (norm(rel_pos(iat,1:2) - pris_pos(iat,1:2)) > 2*max(bondlength))
         for im = -1 : 1
            for jm = -1 : 1
               tmp = pris_pos(iat,:) + im*mcell(1,:) + jm*mcell(2,:);
               if (norm(rel_pos(iat,1:2) - tmp(1:2)) < 2*max(bondlength))
		  pris_pos(iat,:) = tmp;
	       end
            end
         end
      end
   end   
   supercell_pris = zeros(tot_natoms*4,3);
   is = 0;
   for im = [0,-1]
      for jm = [0,-1]
	 for inat = 1 : tot_natoms
	    is = is + 1;
            supercell_pris(is,:) = pris_pos(inat,:) + im*mcell(1,:) + jm*mcell(2,:);
         end
      end
   end
   pris_pos_ws = supercell_pris(IN_ws,:);   
else
   rel_pos_ws = rel_pos;
   pris_pos_ws = rel_pos_ws;
   rel_id_ws = [1:tot_natoms];
end
clear supercell_pris pris_pos

% multilayer or monolayer?
multilayer = false;

if (max(structure.layer)>1)
    multilayer = true;
    outfname = join([outfname,'multilayer_']);
end

twisted = false;
if(any(theta) ~= 0.0)
   twisted = true;
   outfname = join([outfname,'twisted_']);
end
if(multilayer)
   for il = 1 : nlayer-1
      for jl = il + 1 : nlayer
          if(tmdc(il)==tmdc(jl) && bondlength(jl)~=bondlength(il))
             fprintf('WARNING: You have selected the same TMD for layer %i and %i but used different lattice vectors',il,jl)
          end
      end
   end
end

if(task == 5 && bse_irh > sum(natoms))
   error('Error in input: bse_irh index is greater than the number of atoms. Aborting')
end
 
% if SOC double the workspace and change the filename 
nspin = 1;
if (lsoc)
    nspin = 2;
    outfname = join([outfname,'wSOC_']);
end

for ilayer = 2 : nlayer
   outfname = join([outfname,'theta=',num2str(theta(ilayer))]);
end
dirname = join([outfname,'.save']);
if(interlayer_int && Interpd)
   outfname = join([outfname,'_pzdz2']);
end

if(exist(dirname ,'dir'))
   if(~restart)
      fprintf('WARNING: Directory %s exists already. Files will be overwritten! \n\n',dirname)
   else
      % Check if orbital_class.mat exists for restart calculations
      if(~exist(join([dirname,'/orbital_class.mat'])))
         error('Cannot find orbital_class.mat. Aborting restart ...')
      end

      % Check if H.mat exists for restart calculations
      if(~exist(join([dirname,'/H.mat'])) && read_ham)
         error('Cannot find H.mat. Aborting restart ...')
      end

      if(~exist(join([dirname,'/gradH.mat'])) && task==5)
         error('Cannot find gradH.mat. Aborting restart ...')
      end

   end
else
   if(restart)
      error('Cannot find .save folder. Aborting restart ...')
   else
      [status, msg, msgID] = mkdir(dirname);
      if(~status)
         error('Error in main.m: error in creating directory. Aborting ...')
      end
   end
end

% Exit here if only the directory is needed
if(onlydir)
   exit
end

fprintf('--> Setting all irreducible parameters ... ')

if(multilayer)
   if(nlayer ~= max(structure.layer))
       error('Number of layers in position file does not correspond to input')
   end
end
% Set all irreducible parameters
onsite = zeros(11,nlayer);
type1 =  zeros(11,11,nlayer);
type2 =  zeros(11,11,nlayer);
type3 =  zeros(11,11,nlayer);
type4 =  zeros(11,11,nlayer);
type5 =  zeros(11,11,nlayer);
type6 =  zeros(11,11,nlayer);
lsoM = zeros(nlayer,1);
lsoX = zeros(nlayer,1);
for i = 1 : nlayer
   [onsite(:,i), type1(:,:,i), type2(:,:,i), type3(:,:,i), type4(:,:,i), type5(:,:,i), type6(:,:,i),lsoM(i),lsoX(i)] = set_parameters(tmdc(i),gw_par);
end

% Set interlayer parameters
[pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm,pd_vint_lay2_z,pd_vint_lay2_parm]=inter_par(multilayer,nlayer,tmdc);
fprintf('done\n\n')

% Number of orbitals, number of occupied states and number of eigenvalues
for ilayer = 1 : nlayer
   norbs(ilayer) = nspin*natoms(ilayer)/3*11;
   nspinorbs(ilayer) = sum(natoms)/3*11;
end
tot_norbs = sum(norbs);
noccs = nspin*7*tot_natoms/3;

if(reduced_workspace)
    neigs = 2*num_cond; %nspin*7*tot_natoms/3 + num_cond;
else
    neigs = tot_norbs;
end
% Check number of valence states + conduction states is < total number of states
if(neigs > tot_norbs)
   fprintf('Number of conduction states too large, reduce num_cond. max(num_cond) = %i, you have set num_cond = %i \n',tot_norbs-nspin*7*tot_natoms/3,num_cond)
   error('Aborting ...')
end

% if(multilayer)
%    if(twisted)
%       if(tmdc(1)==tmdc(2))
%          if(lsoc)
%             fprintf(['* SOC calculation on ',TMD1,' ', TMD2,' twisted homobilayer *','\n\n'])
%          else
%             fprintf(['* Calculation on ',TMD1,' ', TMD2,' twisted homobilayer *','\n\n'])
%          end
%       else
%          if(lsoc)
%             fprintf(['* SOC calculation on ',TMD1,' ', TMD2,' twisted heterobilayer *','\n\n'])
%          else
%             fprintf(['* Calculation on ',TMD1, ' ', TMD2, ' twisted heterobilayer *','\n\n'])
%          end
%       end
%    else
%       if(tmdc(1)==tmdc(2))
%          if(lsoc)
%             fprintf(['* SOC calculation on ',TMD1,' ', TMD2,' homobilayer *','\n\n'])
%          else
%             fprintf(['* Calculation on ',TMD1, ' ', TMD2, ' homobilayer *','\n\n'])
%          end
%       else
%         if(lsoc)
%             fprintf(['* SOC calculation on ',TMD1,' ', TMD2,' heterobilayer *','\n\n'])
%         else
%             fprintf(['* Calculation on ',TMD1, ' ', TMD2, ' heterobilayer *','\n\n'])
%         end
%       end
%    end
% else
%    if(lsoc)
%       fprintf(['* SOC calculation on ',TMD1, ' ', ' monolayer *','\n\n'])
%    else
%       fprintf(['* Calculation on ',TMD1, ' ', ' monolayer *','\n\n'])
%    end
% end

theta = theta*pi/180.0;
Rtheta = cell(nlayer,1);
for ilayer = 1 : nlayer
   Rtheta{ilayer} = [cos(theta(ilayer)) -sin(theta(ilayer)) 0; sin(theta(ilayer)) cos(theta(ilayer)) 0; 0 0 1];
end

% Size of supercell in real space
if(tot_natoms==3 || tot_natoms==6 || tot_natoms == 9 || tot_natoms == 42)
    ncell1 = 5;
    ncell2 = 5;
else
    ncell1 = 3;
    ncell2 = 3;
end

% Lattice vectors of unit cell
a1 = unit_cell(1,:);
a2 = unit_cell(2,:);
a3 = unit_cell(3,:);
% Unit lattice vector
u1 = unit_cell(1,:)/norm(unit_cell(1,:));
u2 = unit_cell(2,:)/norm(unit_cell(2,:));
u3 = unit_cell(3,:)/norm(unit_cell(3,:));
% Lattice vectors of moire cell or supercell
ma1 = mcell(1,:);
ma2 = mcell(2,:);
ma3 = mcell(3,:);

if(convention == 2)
    a2= a2-a1;
    u2= u2-u1;
    ma2= ma2-ma1;
end

% Reciprocal lattice vectors of unit cell
v=abs(dot(a1,cross(a2,a3)));
b1=2*pi*cross(a2,a3)/v;
b2=2*pi*cross(a3,a1)/v;
b3=2*pi*cross(a1,a2)/v;
% Reciprocal lattice vectors of moire cell
mv=abs(dot(ma1,cross(ma2,ma3)));
mb1=2*pi*cross(ma2,ma3)/mv;
mb2=2*pi*cross(ma3,ma1)/mv;
mb3=2*pi*cross(ma1,ma2)/mv;

if(read_kpts)
   if(exist('kpts.in','file'))
      data = importdata('kpts.in',' ');
      recL=[mb1;mb2;mb3];
      all_kpts = data;
      knum_tot = size(all_kpts,1)
      scale_axis = zeros(knum_tot,1);
      for ik = 2 : knum_tot
          dk = norm(all_kpts(ik,:)-all_kpts(ik-1,:));
          scale_axis(ik) = scale_axis(ik-1) + dk;
      end
   else
     error('Error main.m cannot find kpts.in. Aborting ...')
   end
else
   % Generate k-points
   [all_kpts,scale_axis,knum_tot,recL] = generate_kpoints(task,multilayer,twisted,miniBZ, ...
            mb1,mb2,mb3,b1,b2,b3,knum,bse_shifted);
end

% Check number of workers is not greater than number of k-points
if(num_workers > knum_tot)
   fprintf('num_workers grater than total number of k-points (%i). Please reduce number of k-points', knum_tot)
   error('Aborting ...')
end

T = 0.0;
Hsoc = cell(nlayer,1);
Hmat = cell(1,knum_tot);
if(task == 5)
   gradH_x = cell(1,knum_tot);
end
if(restart)
   fprintf('--> Reading orbitals from orbital_class.mat ... ')
   load(join([dirname,'/orbital_class.mat']));
   % Check number of orbitals is correct
   if(size(orbitals,2) ~= tot_norbs)
      error('Number of orbitals from restart does not match. Aborting ... ')
   end
   fprintf('done\n\n')
   fprintf('--> Reading rotation matrix from unsym_mat.mat ... ')
   load(join([dirname,'/unsym_mat.mat']));
   % Check size of rotation matrix
   if(size(T,2) ~= tot_norbs)
      error('Size of rotation matrix from restart is not correct. Aborting ... ')
   end
   fprintf('done\n\n')
   if(read_ham)
      fprintf('--> Reading Hamiltonian from H.mat ... ')
      load(join([dirname,'/H.mat']));
      fprintf('done\n\n')
      % Check size is correct
      if(size(Hmat{1},1) ~= tot_norbs)
         error('Error in main.m: Size of Hamiltonian matrix is wrong. Aborting ... ')
      end
      if(size(Hmat,2) ~= knum_tot)
         error('Error in main.m: Number of Hamiltonian matrces is not equal to number of k-points. Aborting ... ')
      end
      % If save workspace check if H is sparse
      if(reduced_workspace)
         if(~issparse(Hmat{1}))
            fprintf(' * Hamiltonian matrices are full. Making them sparse ...')
            for ik = 1 : knum_tot
               Hmat{ik} = sparse(Hmat{ik});
            end
            fprintf('done *\n\n')
         end
      else 
         if(issparse(Hmat{1}))
            fprintf(' * Hamiltonian matrices are sparse. Making them full ...')
            for ik = 1 : knum_tot
               Hmat{ik} = full(Hmat{ik});
            end
            fprintf('done *\n\n')
         end
      end
      if(task==5)
         fprintf('--> Reading Hamiltonian gradient from gradH.mat ... ')
         load(join([dirname,'/gradH.mat']));
         fprintf('done\n\n')
         if(size(gradH_x{1},1) ~= tot_norbs)
            error('Error in main.m: Size of Hamiltonian gradient matrix is wrong. Aborting ... ')
         end
         if(size(gradH_x,2) ~= knum_tot)
            error('Error in main.m: Number of Hamiltonian gradient matrces is not equal to number of k-points. Aborting ... ')
         end
         % If save workspace check if H is sparse
         if(reduced_workspace)
            if(~issparse(gradH_x{1}))
               fprintf(' * Hamiltonian gradient matrices are full. Making them sparse ...')
               for ik = 1 : knum_tot
                  gradH_x{ik} = sparse(gradH_x{ik});
               end
               fprintf('done *\n\n')
            end
         else 
            if(issparse(gradH_x{1}))
               fprintf(' * Hamiltonian gradient matrices are sparse. Making them full ...')
               for ik = 1 : knum_tot
                  gradH_x{ik} = full(gradH_x{ik});
               end
               fprintf('done *\n\n')
            end
         end
      end 
   end
else
   fprintf('--> Initialising orbitals ... ')
   % Create an array of orbital classes and initialises it
   orbitals(tot_norbs) = orbital();
   for iorb = 1 : tot_norbs
       orbitals(iorb) = orbital(); 
       orbitals(iorb).Ham_index = 0;
       orbitals(iorb).Found_centre = false;
   end

   % Initialise orbitals
   [orbitals] = initialise_orbitals2(orbitals,structure,rel_pos,...
        rel_pos_ws,pris_pos_ws,rel_id_ws,mcell,tot_natoms,natoms,multilayer,nlayer,nspin);
   fprintf('done\n\n')
   %clear rel_pos pris_pos rel_pos_ws pris_pos_ws rel_id_ws

   fprintf('--> Finding neighbors and assigning hopping parameters ... ')
   % Find neighbors and assign hopping parameters
   [orbitals] = opt_find_neigh2(orbitals,mcell,ncell1,ncell2,bondlength,aXX2,...
      nspin,type1,type2,type3,type4,type5,type6,u1,u2,Rtheta,multilayer,nlayer,interlayer_int,Interpd,flipped);
   fprintf('done\n\n')

   % Find transformation matrix from symmetrised basis to unsymmetrised basis
   if(lsoc || interlayer_int)
      clear T;
      fprintf('--> Generating rotation matrix (unsymmetrised basis) ... ')
      T = generate_T(orbitals,tot_norbs,nspin,nlayer);
      fprintf('done\n\n')
      if(reduced_workspace)
         T = sparse(T);
      end
   end
   if(save_workspace)
      fprintf('--> Saving orbitals class into orbital_class.mat ... ')
      save(fullfile(dirname,'orbital_class'),'orbitals');
      fprintf('done\n\n')
      if(ismatrix(T))
         fprintf('--> Saving rotation matrix into unsym_mat.mat ... ')
         save(fullfile(dirname,'unsym_mat'),'T');
         fprintf('done\n\n')
      end
      % if(Hsoc~=0.0)
      %    fprintf('--> Saving SOC into Hsoc.mat ... ')
      %    save('Hsoc','Hsoc');
      %    fprintf('done\n\n')
      % end
   end
end

if(~read_ham)
   % Compute on-site SOC <phi_{i,\sigma} | L.S | phi_{j,\sigma'}>, i.e. the matrix element
   % of L.S in the Wannier basis. This is an on-site term and therefore does
   % not carry momentum, just needs to be computed once.
   if (lsoc)
       clear Hsoc;
       fprintf('--> Computing SOC ... ')
       for ilayer = 1 : nlayer 
           [Hsoc{ilayer},CGp,lang1] = soc(lsoM(ilayer),lsoX(ilayer));
       end
       fprintf('done\n\n')
       if(reduced_workspace)
          for ilayer = 1 : nlayer
             Hsoc{ilayer} = sparse(Hsoc{ilayer});
          end
       end
   end
   
end

% Start parallel Hamiltonian diagonalization
[Hmat,tb_bands,tb_vecs,gradH_x] = opt_build_and_diag_H(task,Hmat,orbitals,norbs,tot_norbs,...
              multilayer,nspin,nlayer,neigs,...
              all_kpts,knum_tot,interlayer_int,lsoc,T,onsite,compute_eigvecs,...
              pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
              pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,flipped,...
              read_ham,save_workspace,reduced_workspace,num_workers,gradient,lambda,...
              pc,dirname,parallel,write_ham,ef_strength);

%if(task==1 || task==3) 
%   clear orbitals;
%end

% Assume conduction bands energies are all above 0 eV
noccs_k = zeros(knum_tot,1);
if(multilayer && twisted)
   if(reduced_workspace)
      for ik = 1 : knum_tot
         ind = find(tb_bands(:,ik) < 0);
         [m,noccs_k(ik)] = max(tb_bands(ind,ik));
      end
   end
   clear ind
else
   noccs_k(:) = noccs;
end

% Plot bandstructure
if(task==1)
   if(save_workspace)
      fprintf('--> Saving eigenvalues and k-path in %s ... ',outfname)
      save(fullfile(dirname,outfname),'scale_axis','tb_bands');
      fprintf('done\n\n')
   end
   file_ID = fopen(join([outfname,'_BS.dat']),'w');
   for in = 1 : neigs
      for ik = 1 : knum_tot
         fprintf(file_ID,'%1.8f %5.8f\n',scale_axis(ik),tb_bands(in,ik));
      end
      fprintf(file_ID,'\n');
   end
   fclose(file_ID);
   if(~read_kpts)
      plot_bandstructure;
   end
   %clear orbitals;
end

if(task==2)
   phase1 = complex(zeros(tot_norbs,1));
   phase2 = complex(zeros(tot_norbs,1));
   for iorb = 1 : tot_norbs
      if(orbitals(iorb).l==1)
         x = orbitals(iorb).Pcentre;
      elseif(orbitals(iorb).l==2)
         x= orbitals(iorb).Centre;
      end
      phase1(iorb) = exp(-1i*dot(mb1,x));
      phase2(iorb) = exp(-1i*dot(mb2,x));
   end
   clear x
   %clear orbitals
   fprintf('--> Plotting Berry curvature and computing Chern number ... ')
   c = chern(tb_vecs,knum_tot,phase1,phase2,multilayer,twisted,nspin,tot_norbs,noccs_k,all_kpts,mb1,mb2,outfname);
   fprintf('done\n\n')
   fprintf('Chern number of selected bands is %i.\n\n', c)
   clear orbitals;
end

% Plot DOS
if (task==3)
    plot_DOS;
end

% Plot density
if(task==4)
   nr1A=zeros(1,natoms(1));
   nr1B=zeros(1,natoms(1));
   nr2A=zeros(1,natoms(2));
   nr2B=zeros(1,natoms(2));
   fprintf('--> Plotting charge density ... ')
   [nr1A,nr1B] = plot_den(tb_vecs,orbitals,structure.x(1:natoms(1)),structure.y(1:natoms(1)),structure.z(1:natoms(1)), ...
                          structure.name(1:natoms(1)),mcell,natoms(1),tot_norbs,noccs_k,knum_tot,knum,outfname, ...
                          false,true,1,metal1,chalc1,nspin,twisted,multilayer);
   [nr2A,nr2B] = plot_den(tb_vecs,orbitals,structure.x(natoms(1)+1:tot_natoms),structure.y(natoms(1)+1:tot_natoms), ...
                          structure.z(natoms(1)+1:tot_natoms),structure.name(natoms(1)+1:tot_natoms),mcell,natoms(2),...
                          tot_norbs,noccs_k,knum_tot,knum,outfname,false,true,2,metal2,chalc2,nspin,twisted,multilayer);
   fprintf('done\n\n')
   clear nr1A nr1B nr2A nr2B orbitals;
end

if(task==5)
   c_index = cell(tot_natoms,1);
   c_index_up = cell(tot_natoms,1);
   c_index_down = cell(tot_natoms,1);
   % Find indeces of orbitals with same centre
   for inat = 1 : tot_natoms
       centre = [structure.x(inat),structure.y(inat),structure.z(inat)];
       if(structure.name(inat) == "Mo" || structure.name(inat) == "W")
          orb = findobj(orbitals,'Centre',centre);
          c_index{inat} = [orb.Ham_index];
          if(lsoc)
             orb = findobj(orbitals,'Centre',centre,'Spin',1);
             c_index_up{inat} = [orb.Ham_index];
             orb = findobj(orbitals,'Centre',centre,'Spin',-1);
             c_index_down{inat} = [orb.Ham_index];
          else
             c_index_up{inat} = c_index{inat};
             c_index_down{inat} = c_index{inat};
          end
       else
          orb = findobj(orbitals,'Pcentre',centre);
          c_index{inat} = [orb.Ham_index];       
          if(lsoc)
             orb = findobj(orbitals,'Pcentre',centre,'Spin',1);
             c_index_up{inat} = [orb.Ham_index];       
             orb = findobj(orbitals,'Pcentre',centre,'Spin',-1);
             c_index_down{inat} = [orb.Ham_index];
          else 
             c_index_up{inat} = c_index{inat};
             c_index_down{inat} = c_index{inat};
          end
       end
   end
   Cv = complex(zeros(tot_norbs,vn,knum_tot));
   Cc = complex(zeros(tot_norbs,cn,knum_tot));
   Ev = zeros(vn,knum_tot);
   Ec = zeros(cn,knum_tot);
   v = zeros(vn,knum_tot);
   c = zeros(cn,knum_tot);
   % Find indeces of orbitals associated to hole centre
   r_h = [structure.x(bse_irh),structure.y(bse_irh),structure.z(bse_irh)];
   orb_h_up = findobj(orbitals,'Centre',r_h,'Spin',1);
   orb_h_down = findobj(orbitals,'Centre',r_h,'Spin',-1);
   if(isempty(orb_h_up))
      orb_h_up = findobj(orbitals,'Pcentre',r_h,'Spin',1);
      orb_h_down = findobj(orbitals,'Pcentre',r_h,'Spin',-1);
   end
   % Move hole position to the centre of supercell for plotting
   %if(~twisted)
      r_h = r_h + (knum-1)/2*(mcell(1,:) + mcell(2,:));
   %end
   % Valence states
   for ik = 1 : knum_tot
      for iv = 1 : vn
         v(iv,ik) = noccs_k(ik) - (vn-iv);
         Cv(:,iv,ik) = tb_vecs(:,v(iv,ik),ik);
         Ev(iv,ik) = tb_bands(v(iv,ik),ik);
      end
      % Conduction states
      for ic = 1 : cn
         c(ic,ik) = noccs_k(ik) + ic;
         Cc(:,ic,ik) = tb_vecs(:,c(ic,ik),ik);
         Ec(ic,ik) = tb_bands(c(ic,ik),ik);
      end
   end

   clear Hmat orbitals tb_bands tb_vecs

   delete(gcp('nocreate'))
   if(bse_serial)
      [D1,BSE_eig,omega,Resigma,Acvk,osc_strength,ind1,ind2,ind3] =  BSE(vn,cn,Cv,Cc,Ev,Ec,knum_tot,all_kpts,...
                     bondlength(1)*sqrt(3),gradH_x,outfname,mcell,c_index,...
                     tot_natoms,structure,elho_int,dL);
   else
      [BSE_eig,omega,Resigma,Acvk,osc_strength,ind1,ind2,ind3] =  BSE_parallel(vn,cn,Cv,Cc,Ev,Ec,knum_tot,all_kpts,...
                     bondlength(1)*sqrt(3),gradH_x,outfname,mcell,c_index,...
                     tot_natoms,structure,elho_int,dL);
   end
   fprintf('--> Plotting excitonic wavefunction, for exciton %i ...\n', bse_eig_plot)
   [BSE_eig, ind_eig] = sort(BSE_eig);
   Acvk = Acvk(:,ind_eig);
   format long
   BSE_eig(1:20)
   [supercell, psi_M1, psi_M2, psi_kM1, psi_kM2] = plot_exciton_wfc2(c_index_up,c_index_down,r_h, orb_h_up, orb_h_down, Acvk(:,bse_eig_plot:bse_eig_plot+1), Cv, Cc, mcell, ...
                        structure.x, structure.y, structure.z, all_kpts, knum_tot, tot_natoms, ...
                        ind1, ind2, ind3, size(ind1,1),BSE_eig(bse_eig_plot),join(['exciton_wfc_eig',num2str(bse_eig_plot)]), twisted);
   fprintf('done.\n')
   if(save_workspace)
      BSE_outfname = join([outfname,'_BSE']);
      fprintf('--> Saving BSE eigenvalues and oscillator strength in %s ... ',BSE_outfname)
      save(fullfile(dirname,BSE_outfname),'BSE_eig','Resigma','omega','osc_strength', ...
      'Acvk', 'ind1', 'ind2', 'ind3', 'c_index_up', 'c_index_down', 'Cv', 'Cc', ... 
      'psi_M1', 'supercell','psi_kM1','psi_kM2');
      fprintf('done\n\n')
   end
   clear Acvk BSE_eig Resigma omega osc_strength ind1 ind2 ind3 c_index_up c_index_down Cv Cc ...
         psi_M1 psi_kM1 psi_kM2
   delete(gcp('nocreate'))
end

if(task==6)
   fprintf('--> Plotting wavefunction at Gamma, K and M ')
      psink_G = plot_wfc_new(tb_vecs(:,noccs_k(1),1),orbitals,structure.x,...
                structure.y,structure.z,structure.name,mcell,...
                natoms(1)+natoms(2),size(orbitals,2),outfname,twisted,...
                multilayer,nlayer,'Gamma');
      psink_M = plot_wfc_new(tb_vecs(:,noccs_k(knum),knum),orbitals,structure.x,...
                structure.y,structure.z,structure.name,mcell,...
                natoms(1)+natoms(2),size(orbitals,2),outfname,twisted,...
                multilayer,nlayer,'M');
      psink_K = plot_wfc_new(tb_vecs(:,noccs_k(2*knum),2*knum),orbitals,structure.x,...
                structure.y,structure.z,structure.name,mcell,...
                natoms(1)+natoms(2),size(orbitals,2),outfname,twisted,...
                multilayer,nlayer,'K');
   fprintf('done\n\n')
   clear orbitals
end

if(task==7 && write_ham)
   write_ham_hdf5
end

%close all;
%clear tb_vecs tb_bands; 
t = clock;
fprintf('%i %s %i %02i:%02i\n',t(3),months{t(2)},t(1),t(4),t(5))
clear t;
fprintf('*** End calculation ***\n')
toc
