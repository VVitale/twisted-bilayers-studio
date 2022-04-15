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
%%   D1, D2 : Eigenvalue of single-particle and BSE eq.    %%
%%   omega : array with h\omega/2\pi for plotting          %%
%%   Resigma : Real part of sigma_xx                       %%
%%   Acvk : Eigenvectors of BSE eq.                        %%
%%   s : Oscillator strenght                               %%
%%   ind1, ind2, ind3 : Indices for c,v,k in BSE Hamilton. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                         %%
%% Written by Valerio Vitale                               %%
%%                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D1,D2,omega,Resigma,Acvk,s,ind1,ind2,ind3] = BSE(vn,cn,Cv,Cc,Ev,Ec,knum_tot,all_kpts,...
                                                 a,gradH_x,outfname,mcell,index,nat,...
                                                 structure,elho_int,dL)

% Dimension of BSE Hamiltonian 
dim = vn*cn*knum_tot;
% Number of BSE eigenvalues to compute (for the moment compute all of them)
neig = vn*cn*knum_tot;
% Volume of unit cell
V = sqrt(3)*a^2/2*knum_tot;
hbar = 6.582119569*10^(-16);

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

ind1 = ind(:,1);
ind2 = ind(:,2);
ind3 = ind(:,3);

% Compute Ec-Ev for single-particle (no interaction)
s = zeros(dim,1);
for idim = 1 : dim
   D1(idim) = (Ec(ind1(idim),ind3(idim))-Ev(ind2(idim),ind3(idim)));
   s(idim) = abs(Cv(:,ind2(idim),ind3(idim))'*gradH_x{ind3(idim)}*Cc(:,ind1(idim),ind3(idim)))^2;
end

% Compute single-particle absorbance spectrum
D1 = sort(D1);
omega = linspace(0.5,3.5,1000);
Resigma = zeros(length(omega),1);
L = zeros(dim,1);
idim = 0;
for iom = 1 : length(omega)
   for idim = 1 : dim
       L(idim) = lorentzian(omega(iom),D1(idim),dL);
   end
   Resigma(iom) = 4*pi/(omega(iom)*V)*dot(s,L);
end

% Plot single-particle spectrum
plot(omega,Resigma)
clear L

% Save single-particle spectrum
file_ID = fopen(join([outfname,'_OpticalCond_nointer.dat']),'w');
for iom = 1 : length(omega)
    fprintf(file_ID,'%2.8f %4.8f\n',omega(iom),Resigma(iom));
end
fclose(file_ID);

% Compute absorbance spectrum with interaction
D2 = zeros(size(D1));
if(elho_int)
   hold on
   clear D Resigma

   % Compute  Icc and Ivv Eq. 7 in Ref. 1
   Icc = cell(cn*knum_tot,cn*knum_tot);
   Ivv = cell(vn*knum_tot,vn*knum_tot);
   t = zeros(nat,1);
   idim = 0;
   jdim = 0;
   for ic = 1 : cn
      for ik = 1 : knum_tot
         idim = idim + 1;
         jdim = 0;
         for icp = 1 : cn
            for ikp = 1 : knum_tot
               jdim = jdim + 1;
               for inat = 1 : nat
                  t(inat)  = dot(Cc(index{inat},ic,ik),Cc(index{inat},icp,ikp));
               end
               Icc{idim,jdim} = t;
            end
         end
      end
   end
   
   idim = 0;
   jdim = 0;
   t = zeros(nat,1);
   for iv = 1 : vn
      for ik = 1 : knum_tot
         idim = idim + 1;
         jdim = 0;
         for ivp = 1 : vn
             for ikp = 1 : knum_tot
                jdim = jdim + 1;
                for inat = 1 : nat
                    t(inat) = dot(Cv(index{inat},iv,ik),Cv(index{inat},ivp,ikp));
                end
                Ivv{idim,jdim} = t;
            end
         end
      end
   end
   
   % Compute Keldysh potential in real-space
   [VR,lattice] = real_pot2(mcell,structure,a,nat,sqrt(knum_tot));
   
   % Compute Keldysh potential in reciprocal space through lattice Fourier transform
   u = cell(knum_tot,knum_tot); 
   tic
   for ik = 1 : knum_tot
      q = all_kpts(ik,:)-all_kpts;
      for ikp = 1 : knum_tot
          ukkp = screenC2(q(ikp,:),knum_tot,nat,VR,lattice,structure,mcell);
          if(norm(ukkp - ukkp')<1E-8)
             u{ik,ikp} = ukkp;
          else
             norm(ukkp - ukkp')
             error('The potential is not Hermitian')
          end
      end
   end
   toc
   
   tic
   % Build BSE Hamiltonian
   W = zeros(dim);
   for jdim = 1:dim
       icp = ind1(jdim);
       ivp = ind2(jdim);
       ikp = ind3(jdim);
       icp_loc = (icp-1)*knum_tot + ikp;
       ivp_loc = (ivp-1)*knum_tot + ikp;
       v = zeros(dim,1);
       for idim = 1 : dim
           ic = ind1(idim);
           iv = ind2(idim);
           ik = ind3(idim);
           ic_loc = (ic-1)*knum_tot + ik;
           iv_loc = (iv-1)*knum_tot + ik;
           if(idim==jdim)
              v(idim) = v(idim) + Ec(ic,ik) - Ev(iv,ik);
           end
           v(idim) = v(idim) + 1/knum_tot*(Icc{icp_loc,ic_loc}'*(u{ik,ikp}*Ivv{ivp_loc,iv_loc}));
      end
      W(:,jdim) = v(:);
   end
   toc
   clear Icc Ivv
   
   % Check W is symmetric
   if(norm(W' - W) > 1E-8)
      fprintf('hole-electron matrix is not Hermitian')
   end

   % Diagonalise BSE Hamiltonian
   [Acvk,D2] = eig(W);
   D2 = diag(D2);
   
   % In neig != dim pad 
   Acvk = [Acvk,eye(dim,dim-neig)];
   D2 = [D2;1000*ones(dim-neig,1)];
   

   % Compute absorbance spectrum 
   s = zeros(neig,1);
   sumM = 0;
   for idim = 1 : neig
      sumM =0;
      for jdim = 1 : dim
          sumM = sumM + Acvk(jdim,idim)*(Cv(:,ind2(jdim),ind3(jdim))'*gradH_x{ind3(jdim)}*Cc(:,ind1(jdim),ind3(jdim)));
      end
      s(idim) = abs(sumM)^2;
   end
   for iom = 1 : length(omega)
      L = zeros(neig,1);
      for idim = 1 : neig
         L(idim) = lorentzian(omega(iom),D2(idim),dL);
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
   rho1 = zeros(size(E,2),1);
   rho2 = zeros(size(E,2),1);
   for ie = 1 : size(rho1)
      for idim = 1 : dim
         rho1(ie) = rho1(ie) + lorentzian(E(ie),D1(idim),dL);
         rho2(ie) = rho2(ie) + lorentzian(E(ie),D2(idim),dL);

      end
   end
   
   rho1 = rho1/cn/vn/knum_tot;
   rho2 = rho2/cn/vn/knum_tot;
   % Plot Joint density of states
   figure
   plot(E,rho1,'b',E,rho2,'r')
   saveas(gca,join([outfname,'_JDOS.png']))
   file_ID = fopen(join([outfname,'_JDOS_with_inter.dat']),'w');
   for iom = 1 : length(E)
       fprintf(file_ID,'%2.8f %4.8f\n',E(iom),rho2(iom));
   end

   fclose(file_ID);
   delete(gcp('nocreate'))
   %clear Acvk Cc Cv 
   clear W rho E L rho1 rho2
   clear Ec Ev
end

end
