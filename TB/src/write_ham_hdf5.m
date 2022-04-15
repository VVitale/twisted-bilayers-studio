%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Routine to write Hamiltonian in hdf5      %%
%% format. It also writes info on the        %%
%% moire geometry and k vectors.             %%
%%                                           %%
%% The hamiltonian is divided into real anc  %%
%% imaginary part as matlab hdf5 cannot      %%
%% deal with complex numbers                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                           %%
%% Written by Valerio Vitale, April 2021     %%
%%                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hamiltonian is written in the following format:
% leading dimension (ncolumns) is total number of orbitals * number of k-points
% second dimension (nrows) is total number of orbitals (Norb)
% H =
%      __Norb__________________  ________
%      |      |      |      |     |      |
% Norb{| H(k1)| H(k2)| H(k3)| ... | H(kN)|
%      |      |      |      |     |      |
%      |______|______|______|__  _|______|


% First transform sparse Hamiltonian into full matrix
Hmat_full = cell(knum_tot);
for ik = 1 : knum_tot
   Hmat_full{ik} = full(Hmat{ik});
end

% Create an entry for real and complex part of Hamiltonian
h5create(ham_fname,'/hamiltonian/realH',[tot_norbs,tot_norbs*knum_tot])
h5create(ham_fname,'/hamiltonian/cmplxH',[tot_norbs,tot_norbs*knum_tot])
count = [tot_norbs,tot_norbs];
for ik = 1 : knum_tot
   start = [1,tot_norbs*(ik-1)+1];
   h5write(ham_fname,'/hamiltonian/realH',real(Hmat_full{ik}),start,count)
   h5write(ham_fname,'/hamiltonian/cmplxH',imag(Hmat_full{ik}),start,count)
end

% Create a geometry/crystal entry
h5create(ham_fname,'/crystal/moire_cell',size(mcell))
h5write (ham_fname,'/crystal/moire_cell',mcell)
h5create(ham_fname,'/crystal/n_atoms',1)
h5write (ham_fname,'/crystal/n_atoms',tot_natoms)
h5create(ham_fname,'/crystal/n_orbs',1)
h5write (ham_fname,'/crystal/n_orbs',tot_norbs)
h5create(ham_fname,'/crystal/positions',[tot_natoms 3])
h5write (ham_fname,'/crystal/positions',[structure.x(:) structure.y(:) structure.z(:)])
h5create(ham_fname,'/crystal/kgrid',[knum_tot 3])
h5write (ham_fname,'/crystal/kgrid',all_kpts)

% Create an entry for orbitals properties
% Index as defined in Fang et al. PRB 92 205108 (2015)
% dxz = 1, dyz = 2, pz(bottom) = 3, px(bottom) = 4, py(bottom) = 5
% dz2 = 6, dxy = 7, dx2-y2 = 8, pz(top) = 9, px(top) = 10, py(top) = 11
% N.B. Same orbitals with different spin projections have same index
% So the relative index runs from 1 to 11
h5create(ham_fname,'/orbitals/pattern',[1 tot_norbs])
orb_pattern = [orbitals.Rel_index];
h5write (ham_fname,'/orbitals/pattern',orb_pattern)

% z-component of Spin for each orbital
h5create(ham_fname,'/orbitals/spin',[1 tot_norbs])
spins = [orbitals.Spin];
h5write (ham_fname,'/orbitals/spin',spins)

% Layer to which each orbital belongs to
h5create(ham_fname,'/orbitals/layer',[1 tot_norbs])
layers = [orbitals.Layer];
h5write (ham_fname,'/orbitals/layer',layers)

clear count start Hmat Hmat_full orb_pattern spins
