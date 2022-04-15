%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to build and                 %%
%% diagonalise the TB Hamiltonian        %%
%% in parallel.                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Written by Valerio Vitale, Sept. 2020 %%
%%                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Hmat,tb_bands,tb_vecs,gradH_x] = opt_build_and_diag_H(task,Hmat,orbitals,norbs,tot_norbs,...
              multilayer,nspin,nlayer,neigs,...
              all_kpts,knum_tot,interlayer_int,lsoc,T,onsite,compute_eigvecs,...
              pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
              pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,flipped,...
              read_ham,save_workspace,reduced_workspace,num_workers,gradient,lambda,...
              pc,dirname,parallel,write_ham,ef_strength)

if(~read_ham)
   diagH = zeros(tot_norbs);
   temp_norbs = 0;
   %Set on-site terms
   for layer = 1 : nlayer
      for iorb = 1:norbs(layer)
          orb1 = orbitals(iorb + temp_norbs);
          diagH(orb1.Ham_index,orb1.Ham_index) = diagH(orb1.Ham_index,orb1.Ham_index) + onsite(orb1.Rel_index,layer);
          % Add EF as position-dependent on-site term
	  if(ef_strength ~= 0.0)
	     if(orb1.l == 2)
	        z_coor = dot(orb1.Centre,[0 0 1]);
	     elseif(orb1.l == 1)
	        z_coor = dot(orb1.Pcentre,[0 0 1]);
	     end
             diagH(orb1.Ham_index,orb1.Ham_index) = diagH(orb1.Ham_index,orb1.Ham_index) + ef_strength*z_coor;
          end
      end
      temp_norbs = sum(norbs(1:layer));
   end
   clear orb1;

   % Generate mask for 2H stackings
   Mxz = generate_mask(orbitals,tot_norbs,nlayer,flipped);
   if(reduced_workspace)
      Mxz = sparse(Mxz);
   end

   % Start parallel section
   if(parallel)
      parpool(pc,num_workers);
   else
      evalc('parpool(pc,num_workers)');
   end

   fprintf('--> Starting parallel diagonalisation of TB Hamiltonian on %i workers ... \n',num_workers)
   
   % Diagonalisation Parallel loop over k-points
   spmd
      N1 = neigs; N2 = tot_norbs; N3 = knum_tot;
      globalSize_bands = [N1 N3];
      codistr_bands = codistributor1d(2, codistributor1d.unsetPartition, globalSize_bands);
      localSize_bands = [N1, codistr_bands.Partition(labindex)];
      tb_bands_loc = zeros(localSize_bands);
      globalInd = codistr_bands.globalIndices(2);
      if(compute_eigvecs)
         globalSize_vecs = [N2 N1 N3];
         codistr_vecs = codistributor1d(3, codistributor1d.unsetPartition, globalSize_vecs);
         localSize_vecs = [N2, N1, codistr_vecs.Partition(labindex)];
         tb_vecs_loc = zeros(localSize_vecs);
      end
      if(write_ham)
         codistr_cell = codistributor1d(2, codistributor1d.unsetPartition, [1, N3]);
         Hmat_loc = cell(1,codistr_cell.Partition(labindex));
      end
      if(gradient)
         codistr_cell = codistributor1d(2, codistributor1d.unsetPartition, [1, N3]);
         gradH_x_loc = cell(1,codistr_cell.Partition(labindex));
      end
      for loc_ik = 1 : length(globalInd)
        ik = globalInd(loc_ik);
        if(reduced_workspace)
           H = sparse(diagH);
           gH = sparse(zeros(tot_norbs,tot_norbs));
        else
          H = diagH;
          gH = zeros(tot_norbs);
        end
         k = all_kpts(ik,:);
         % Build Hamiltonian matrix
         %H = cell2mat(Hmat_loc{ik});
         if(gradient)
            % = sparse(complex(zeros(size(Hmat_loc(ik)))));
            % Build Hamiltonian matrix and gradient along kx
            [H,gH] = build_H2(H,orbitals,Mxz,k, ...
                                    nlayer,multilayer,nspin,interlayer_int,lsoc,T, ...
                                    pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
                                    pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,gradient,gH); 
         else
            % Build Hamiltonian matrix and gradient along kx
            H = build_H2(H,orbitals,Mxz,k, ...
                                    nlayer,multilayer,nspin,interlayer_int,lsoc,T, ...
                                    pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm, ...
           			 pd_vint_lay2_z,pd_vint_lay2_parm,theta,Hsoc,Interpd,gradient);
         end

         if(write_ham)
            Hmat_loc(:,loc_ik) = mat2cell(H,tot_norbs,tot_norbs);
         end
         if(gradient)
            gradH_x_loc(:,loc_ik) = mat2cell(gH,tot_norbs,tot_norbs);
         end
         
         if(task~=7)
            if(reduced_workspace)
               % Shift eigenspectrum
               H = H + lambda*sparse(eye(size(H)));
               if (compute_eigvecs)
                  [tb_vecs_loc(:,:,loc_ik),D] = eigs(H,neigs,'SM');
                  [tb_bands_loc(:,loc_ik),ind] = sort(diag(real(D)));
                  tb_vecs_loc(:,:,loc_ik) = tb_vecs_loc(:,ind,loc_ik);
               else
                  D = eigs(H,neigs,'SM');
                  [tb_bands_loc(:,loc_ik),ind] = sort(real(D));
               end       
            else
               if (compute_eigvecs)
                  [V,D] = eig(H,'vector');
                  [tb_bands_loc(:,loc_ik),ind] = sort(real(D));
                  tb_vecs_loc(:,:,loc_ik) = V(:,ind);
               else
                  D = eig(H,'vector');
                  [tb_bands_loc(:,loc_ik),ind] = sort(real(D));
               end
            end
         end
         fprintf('k-vector number %i completed \n',ik);
      end
      tb_bands_global = codistributed.build(tb_bands_loc, codistr_bands);
      if(compute_eigvecs)
         tb_vecs_global = codistributed.build(tb_vecs_loc, codistr_vecs);
      end
      if(write_ham)
         Hmat_global = codistributed.build(Hmat_loc, codistr_cell);
      end
      if(gradient)
         gradH_x_global = codistributed.build(gradH_x_loc, codistr_cell);
      end
   end
   clear tb_bands_loc
   if(compute_eigvecs)
      clear tb_vecs_loc
   end
   % END OF PARALLEL SECTION
   fprintf('done\n\n')

   if(write_ham)
      Hmat = gather(Hmat_global);
      clear Hmat_global Hmat_loc;
   else
      Hmat = 0.0;
   end
   if(gradient)
      gradH_x = gather(gradH_x_global);
      clear gradH_x_global gradH_x_loc;
   else
      gradH_x = 0.0;
   end
   if(save_workspace)
      fprintf('--> Saving Hamiltonian matrix in H.mat ... ')
      save(fullfile(dirname,'H'),'Hmat');
      fprintf('done\n\n')
   end
   
   if(save_workspace && task==5)
      fprintf('--> Saving Hamiltonian gradient matrix in gradH.mat ... ')
      save(fullfile(dirname,'gradH'),'gradH_x');
      fprintf('done\n\n')
   end
   clear diagH H gH 

else
   fprintf('--> Starting parallel diagonalisation of TB Hamiltonian on %i workers ... ',num_workers)

% Start parallel section
if(parallel)
   parpool(pc,num_workers);
else
   evalc('parpool(pc,num_workers)');
end

tb_bands_global=zeros(neigs,knum_tot,'distributed');
if(compute_eigvecs)
   tb_vecs_global=zeros(tot_norbs,neigs,knum_tot,'distributed');
end

   gradH_x = 0.0;
   Hmat_loc = distributed(Hmat);
   spmd
      for ik = drange(1:knum_tot)
         H = cell2mat(Hmat_loc(ik));
         if(reduced_workspace)
            H = H + lambda*sparse(eye(size(H)));
            if (compute_eigvecs)
               [tb_vecs_global(:,:,ik),D] = eigs(H,neigs,'SM');
               [tb_bands_global(:,ik),ind] = sort(diag(real(D)));
               tb_vecs_global(:,:,ik) = tb_vecs_global(:,ind,ik);
            else
               D = eigs(H,neigs,'SM');
               [tb_bands_global(:,ik),ind] = sort(real(D));
            end       
         else
            if (compute_eigvecs)
               [V,D] = eig(H,'nobalance','vector');
               [tb_bands_global(:,ik),ind] = sort(real(D));
               tb_vecs_global(:,:,ik) = V(:,ind);
            else
               D = eig(H,'nobalance','vector');
               [tb_bands_global(:,ik),ind] = sort(real(D));
            end
         end
         fprintf('k-vector number %i completed \n',ik)
      end
   end
   % END OF PARALLEL SECTION
   fprintf('done\n\n')
   clear Hmat_loc H 
end

clear D;

if(~reduced_workspace && compute_eigvecs)
    clear V;
end
if(ismatrix(T)) 
   clear T; 
end
if(lsoc && iscell(Hsoc))
   clear Hsoc
end

% GATHER RESULTS
tb_bands = gather(real(tb_bands_global));
clear tb_bands_global;
if(compute_eigvecs)
   tb_vecs = gather(tb_vecs_global);
   clear tb_vecs_loc;
else
   tb_vecs = 0.0;
end

% Suppress output from parpool
%if(task~=5)
if(parallel)
   delete(gcp);
else
   evalc('delete(gcp)');
end
%end
end

