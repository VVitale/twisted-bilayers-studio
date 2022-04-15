%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to construct the BSE                           %%
%% Hamiltonian and diagonalise it                          %%
%% in tight-binding basis given in :                       %%
%% (1) Ridolfi et al. PRB 97, 205409 (2018)                %%
%% and (2) MacDonald et al. PRB 91, 075310                 %%
%% (2015)                                                  %%
%% - INPUT VARIABLES                                       %%
%%   vn, cn : Number of valence and conduction bands       %%
%%   Cv, CC : Coeffients of TB for valence and cond bands  %%
%%   Ev, Ec : Eigenvalues of TB                            %%
%%   knum_tot : total number of k-points in k-grid         %%
%%   all_kpts : Cartesian coordinates of k-points          %%
%%   a : lattice constant                                  %%
%%   gradH_x : Gradient wrt k_x of Hamiltonian matrix      %%
%%   outfname : Name of output file (string)               %%
%%   mcell : Lattice cell vectors                          %%
%%   nat : Total number of atoms                           %%
%%   structure : Cartesian coordinates of all atoms        %%
%%   elho_int : Wheter to include el-ho interaction        %%
%%   dL : broadening of Lorentzian for plotting            %%
%%                                                         %%
%% - OUTPUT                                                %%
%%   BSE_eig, D1 : Eigenvalue of BSE w/o interaction       %%
%%   omega : array with h\omega/2\pi for plotting          %%
%%   Resigma : Real part of sigma_xx                       %%
%%   Acvk : Eigenvectors of BSE                            %%
%%   s : Oscillator strenght                               %%
%%   indc, indv, indk : Indices for c,v,k in BSE Hamilton. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                         %%
%% Written by Valerio Vitale                               %%
%%                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BSE_eig,omega,Resigma,Acvk,s,indc,indv,indk] = BSE_parallel4(vn,cn,Cv,Cc,Ev,Ec,knum_tot,all_kpts,...
                                                 a,gradH_x,outfname,mcell,index,nat,...
                                                 positions,elho_int,dL)

% Dimension of BSE Hamiltonian 
dim = vn*cn*knum_tot;
% Number of BSE eigenvalues to compute (for the moment compute all of them)
neig = 50;% vn*cn*knum_tot;
% Volume of unit cell
V = sqrt(3)*a^2/2*knum_tot;
hbar = 6.582119569*10^(-16);
bse_thr = 1E-10;

% Fix gauge of Cv and Cc by imposing their sum to be
% a real number
for iv = 1 : vn
   for ik = 1 : knum_tot
      s = sum(Cv(:,iv,ik));
      phi = atan2(imag(s),real(s));
      phase = exp(-1i*phi);
      Cv(:,iv,ik) = phase.*Cv(:,iv,ik);
      if(imag(sum(Cv(:,iv,ik)))>10^(-10))
          [iv,ik]
          error('The sum of the valence band coefficients is not real.')
      end
   end
end
for ic = 1 : cn
   for ik = 1 : knum_tot
      s = sum(Cc(:,ic,ik));
      phi = atan2(imag(s),real(s));
      phase = exp(-1i*phi);
      Cc(:,ic,ik) = phase*Cc(:,ic,ik);
      if(imag(sum(Cc(:,ic,ik)))>10^(-10))
          [ic,ik]
          error('The sum of the conduction band coefficients is not real.')
      end
   end
end

% Define indices for conduction and valence bands and k-points
ind = zeros(dim,3);
idim = 0;
for ic = 1 : cn
   for iv = 1 : vn
       for ik = 1 : knum_tot
          idim = idim + 1;
          ind(idim,:) = [ic,iv,ik];
       end
   end
end

indc = ind(:,1);
indv = ind(:,2);
indk = ind(:,3);

omega = linspace(0.5,3.5,1000);
Resigma = zeros(length(omega),1);

% Compute absorbance spectrum with interaction
BSE_eig = zeros(neig,1);
if(elho_int)
   hold on
   clear D Resigma

   fprintf('--> Starting parallel construction of BSE Hamiltonian ...\n')

   % Compute lattice Fourier transform of Keldysh potential
   u0 = screenC4([0,0,0],sqrt(knum_tot),mcell,a);

   % Compute diagonal part
   fprintf('Computing diagonal part of BSE Hamiltonian ...')
   diag_H_BSE = zeros(dim,1);
   parfor idim = 1 : dim
      ic = indc(idim);
      iv = indv(idim);
      ik = indk(idim);
      Icc = zeros(nat,1);
      Ivv = zeros(nat,1);
      for inat = 1 : nat
              Icc(inat) = dot(Cc(index{inat},ic,ik),Cc(index{inat},ic,ik));
              Ivv(inat) = dot(Cv(index{inat},iv,ik),Cv(index{inat},iv,ik));
      end
      diag_H_BSE(idim) = Ec(ic,ik) - Ev(iv,ik) + Icc'*u0*Ivv/knum_tot;
   end
   fprintf('done.\n')

   fprintf('Pre-computing reciprocal-space Keldysh potential ...\n')
   ukkp = complex(zeros(knum_tot));% u0*eye(knum_tot);

   parfor ik = 1 : knum_tot
       for ikp = 1 : knum_tot
          q = all_kpts(ik,:) - all_kpts(ikp,:);
          ukkp(ik,ikp) = screenC4(q,sqrt(knum_tot),mcell,a);
       end
       fprintf('k-point %i of %i done\n',ik,knum_tot)
   end
   fprintf('done.\n')

   tic
   H_BSE = sparse(complex(zeros(dim)));
   parfor jdim = 1 : dim
       icp = indc(jdim);
       ivp = indv(jdim);
       ikp = indk(jdim);
       v = complex(zeros(dim,1));
       for idim = jdim + 1 : dim
          ic = indc(idim);
          iv = indv(idim);
          ik = indk(idim);

          Icc = zeros(nat,1);
          Ivv = zeros(nat,1);

          for inat = 1 : nat
                  Icc(inat) = dot(Cc(index{inat},icp,ikp),Cc(index{inat},ic,ik));
                  Ivv(inat) = dot(Cv(index{inat},ivp,ikp),Cv(index{inat},iv,ik));
          end
          v(idim)= Icc'*ukkp(ik,ikp)*Ivv/knum_tot;
      end
      v(abs(v) < bse_thr ) = complex(0.0,0.0);
      H_BSE(:,jdim) = sparse(v(:));
      fprintf('Row %i of %i completed \n', jdim, dim) 
   end
   toc

   diag(H_BSE)
   H_BSE = H_BSE + H_BSE' + sparse(diag(diag_H_BSE));
   fprintf('done.\n')
   clear Icc Ivv %H_BSE_global
   
   % Check H_BSE is symmetric
   if(norm(full(H_BSE)' - full(H_BSE)) > 1E-8)
      %norm(H_BSE' - H_BSE)
      error('ERROR: hole-electron matrix is not Hermitian. Aborting!')
   end

   % Diagonalise BSE Hamiltonian
   fprintf('--> Starting diagonalization of BSE Hamiltonian ...\n')
   [Acvk,BSE_eig] = eigs(H_BSE,neig,'SR');
   BSE_eig = diag(BSE_eig);
   
   fprintf('done.\n')
   
   % In neig != dim pad 
   Acvk = [Acvk,eye(dim,dim-neig)];
   BSE_eig = [BSE_eig;1000*ones(dim-neig,1)];
   

   % Compute absorbance spectrum 
   s = zeros(neig,1);
   sumM = 0;
   parfor idim = 1 : neig
      sumM =0;
      for jdim = 1 : dim
          sumM = sumM + Acvk(jdim,idim)*(Cv(:,indv(jdim),indk(jdim))'*gradH_x{indk(jdim)}*Cc(:,indc(jdim),indk(jdim)));
      end
      s(idim) = abs(sumM)^2;
   end
   for iom = 1 : length(omega)
      L = zeros(neig,1);
      for idim = 1 : neig
         L(idim) = lorentzian(omega(iom),BSE_eig(idim),dL);
      end
      Resigma(iom) = 4*pi/(V*omega(iom))*dot(s,L);
   end
   
   % Plot absorbance spectrum
   plot(omega,Resigma,'r')

   % Save absorbance spectrum
   saveas(gca,join([outfname,'Re_sigma_xx.png']))
   file_ID = fopen(join([outfname,'_OpticalCond_with_inter.dat']),'w');
   for iom = 1 : length(omega)
       fprintf(file_ID,'%2.8f %4.8f\n',omega(iom),Resigma(iom));
   end
   fclose(file_ID);
   
   % Compute Joint density of states
   E = linspace(0,6,1000);
   rho2 = zeros(size(E,2),1);
   for ie = 1 : size(rho2)
      for idim = 1 : dim
         %rho1(ie) = rho1(ie) + lorentzian(E(ie),D1(idim),dL);
         rho2(ie) = rho2(ie) + lorentzian(E(ie),BSE_eig(idim),dL);

      end
   end
   
   rho2 = rho2/cn/vn/knum_tot;
   % Plot Joint density of states
   figure
   plot(E,rho2,'r')
   saveas(gca,join([outfname,'_JDOS.png']))
   file_ID = fopen(join([outfname,'_JDOS_with_inter.dat']),'w');
   for iom = 1 : length(E)
       fprintf(file_ID,'%2.8f %4.8f\n',E(iom),rho2(iom));
   end

   fclose(file_ID);
   % delete(gcp('nocreate'))
   %clear Acvk Cc Cv 
   clear H_BSE rho E L rho1 rho2
   clear Ec Ev
end

end
