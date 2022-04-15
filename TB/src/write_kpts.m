task = 1
bilayer = true
multilayer = false
knum = 20
path = 'MK' % GM: Gamma-M, MK : M-K, KG: K-Gamma

% Create kpts.in file
fileID = fopen('kpts.in','w');


% Copy-paste unit lattice vectors from the position file into mcell
ucell = [3.18300000 0.00000000 0.00000000
-1.59150000 2.75655886 0.00000000
0.00000000 0.00000000 20.00000000];
ua1 = ucell(1,:);
ua2 = ucell(2,:);
ua3 = ucell(3,:);

v=abs(dot(ua1,cross(ua2,ua3)));
b1=2*pi*cross(ua2,ua3)/v;
b2=2*pi*cross(ua3,ua1)/v;
b3=2*pi*cross(ua1,ua2)/v;

% Copy-paste the moire lattice vectors from the position file into mcell
mcell = [3.18300000 0.00000000 0.00000000
-1.59150000 2.75655886 0.00000000
0.00000000 0.00000000 20.00000000];
ma1 = mcell(1,:);
ma2 = mcell(2,:);
ma3 = mcell(3,:);

% Reciprocal lattice vectors of moire cell
mv=abs(dot(ma1,cross(ma2,ma3)));
mb1=2*pi*cross(ma2,ma3)/mv;
mb2=2*pi*cross(ma3,ma1)/mv;
mb3=2*pi*cross(ma1,ma2)/mv;

[all_kpts,scale_axis,knum_tot,recL,frac_k] = generate_kpoints(task,multilayer,twisted,false,mb1,mb2,mb3,b1,b2,b3,knum);

% Print first 20 k-points
if(strcmp(path,'GM'))
   kinit = 1;
   kfinal = knum;
elseif(strcmp(path,'MK'))
   kinit = knum+1;
   kfinal = 2*knum;
elseif(strcmp(path,'KG'))
   kinit = 2*knum + 1;
   kfinal = 3*knum;
else
   error('k-path not recognised. Aborting')
end

for ik = kinit : kfinal
   fprintf(fileID,'%2.8f  %2.8f  %2.8f\n',all_kpts(ik,:));
end

% Close file
fclose(fileID);
