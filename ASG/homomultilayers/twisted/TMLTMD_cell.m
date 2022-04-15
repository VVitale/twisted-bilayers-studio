%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the moire cell for twisted bilayer graphene (TBLG) given two 
% integers n and m. Here, we assume that m=n+1.
%
% It generates five output files: 
% 	1) XYZ or XSF file for visualization
%	2) .dat file that is used as input file for TB_free.x
%   3) kpoints.in file to be used as input for TB_free.x
%   4) kpaths.dat file to be used for plotting the band structure on the
%      chosen path in k-space
%	5) lammps geometry file
%
% It returns the xyz coordinates of bottom layer (pos) and top layer (pos2) 
% and the total number of atoms (nat).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nat,pos] = TQLTMD_cell(theta_vec,n, a, z, atomlabels, ...
    plot_f, kpath, nk, write_positions, write_cart, write_lammps_input,...
    write_kpath,stackings,orientations,nlayers)
   % Setting parameters 
   m=n+1;
   mass = zeros(3*nlayers,1);
   num_at = zeros(3*nlayers,1);
   for isublayer = 1 : 3*nlayers
       % Define masses and atomic numbers
       if(strcmp(atomlabels{isublayer},'Mo'))
           % Atomic masses
           mass(isublayer) = 95.94;
           % Atoms info
           num_at(isublayer) = 42;
       elseif(strcmp(atomlabels{isublayer},'W'))
           mass(isublayer) = 183.84;
           num_at(isublayer) = 74;
       elseif(strcmp(atomlabels{isublayer},'S1') || strcmp(atomlabels{isublayer},'S2') )
           mass(isublayer) = 32.06;
           num_at(isublayer) = 16;
       elseif(strcmp(atomlabels{isublayer},'Se1') || strcmp(atomlabels{isublayer},'Se2'))
           mass(isublayer) = 78.96;
           num_at(isublayer) = 34;
       end
   end
   mass
   num_at
   
   theta0 = acos((n^2 + 4*n*m + m^2)/(2*(n^2 + n*m + m^2)));
   theta = theta0*theta_vec; 
   
   filename = join(['t',num2str(nlayers),'LTMD']);
   for ilayer = 1 : nlayers
    filename = join([filename atomlabels{(ilayer-1)*3 + 1}],"_");
    filename = join([filename atomlabels{(ilayer-1)*3 + 2} "2"],'');
   end

   % total number of atoms in Moire cell
   tot_nat= 3*nlayers*(n^2 + n*m + m^2);
      
   % Compute Moire supercell  
   u=0.6226;
   c=2*z;
   
   % Lattice vectors of top and bottom layer  
   ualat1 = a*[sqrt(3)/2, sqrt(3)/2;  -0.5,0.5];
   Rphi = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)];
   alat1 = Rphi*ualat1;
   % Lattice vectors of Moire cell
   alatL = [n, m; -m, n+m]*alat1'
   ialatL = inv(alatL');
   
   % Generate reciprocal lattice vectors
   recL = 2*pi*inv(alatL');
   irecL = inv(recL);
   trecL = recL';
   ralat1 = 2*pi*inv(alat1');
   %ralat2 = 2*pi*inv(alat2');

   % Print k-point grid/path to file
%      theta = theta*180/pi
%    str_theta = num2str(theta);
%    str_m = num2str(m);
%    str_n = num2str(n);
%    if (write_kpath)
%       msg = ['Generating kpoint grid for moirè lattice with a twist angle of ',str_theta,' ...'];
%       disp(msg)
%       disp(' ')
%       write_kpts(nk,ralat1,recL,str_theta,str_m,str_n,kpath,1);
%       disp(' ')
%       disp('... done')
%    end
   
   % number of supercell for the search
   period = norm(alatL,1);
   tnat = floor(period/a);
   pos = -10^5*ones(tot_nat,3);
   id(1,tot_nat) = "X";
   delta =  1*10^(-5);
   
%    disp(' ')
%    msg = ['Generating Cartesian atomic coordinates for moirè lattice with a twist angle of ',str_theta,' ...'];
%    disp(msg)
   
   k1=1;
   
   for ilayer = 1 : nlayers
       st=sin(theta(ilayer));
       ct=cos(theta(ilayer));
       Rtheta = [ct -st; st ct];
       for nj = -3*tnat:3*tnat
           for ni = -3*tnat:3*tnat
               r1 = ni*alat1(:,1) + nj*alat1(:,2);
               r2 = r1 + 1/3*alat1(:,1) + 1/3*alat1(:,2);
               %r3 = r2;
               f13 = 1/3;
               f23 = 2/3;
               if(strcmp(orientations{ilayer},'r'))
                   r1 = -r1 + 1/3*(alat1(:,1)+alat1(:,2));
                   r2 = -r2 + 1/3*(alat1(:,1)+alat1(:,2));
                   %r3 = -r2;
                   f13 = 2/3;
                   f23 = 1/3;
               end
               if(strcmp(stackings{ilayer},'t0'))
                   r1r = Rtheta*r1;
                   r2r = Rtheta*r2;
                   r3r = Rtheta*r2;
               elseif(strcmp(stackings{ilayer},'t13'))
                   r1r = Rtheta*(r1+f13*(alat1(:,1)+alat1(:,2)));
                   r2r = Rtheta*(r2+f13*(alat1(:,1)+alat1(:,2)));
                   r3r = r2r;
               elseif(strcmp(stackings{ilayer},'t23'))
                   r1r = Rtheta*(r1+f23*(alat1(:,1)+alat1(:,2)));
                   r2r = Rtheta*(r2+f23*(alat1(:,1)+alat1(:,2)));
                   r3r = r2r;
               end
               l1 = ialatL*r1r;
               l2 = ialatL*r2r;
               l3 = ialatL*r3r;
               if (l1(1) < (1-delta) && l1(1) >= -delta && l1(2) < (1-delta) && l1(2) >= -delta)
                   pos(k1,1)=r1r(1);
                   pos(k1,2)=r1r(2);
                   pos(k1,3)=((ilayer-1)*0.5 + 0.75)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 1};
                   k1 = k1 + 1;
               end
               if (l2(1) < (1-delta) && l2(1) >= -delta && l2(2) < (1-delta) && l2(2) >= -delta)
                   pos(k1,1)=r2r(1);
                   pos(k1,2)=r2r(2);
                   pos(k1,3)= ((ilayer-1)*0.5 + u)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 2};
                   k1 = k1 + 1;
              % end
              % if (l3(1) < (1-delta) && l3(1) >= -delta && l3(2) < (1-delta) && l3(2) >= -delta)
                   pos(k1,1)=r3r(1);
                   pos(k1,2)=r3r(2);
                   pos(k1,3)=((ilayer-1)*0.5+1.5-u)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 3};
                   k1 = k1 + 1;
               end
           end
       end
   end
   nat = k1 - 1;
   pos(:,3) = pos(:,3) - 5;
%Check number of atoms is 3*(n^2 + n*m + m^2)
   if (nat ~= 3*nlayers*(n^2 + n*m + m^2))
       nat
       3*nlayers*(n^2 + n*m + m^2)
       error('Wrong number of atoms');
   end
   disp('... done')
   
   %figure
   theta = theta*180/pi;
   for ilayer = 1 : nlayers
       if(strcmp(stackings(ilayer),'t0') || strcmp(stackings(ilayer),'t13') || ...
               strcmp(stackings(ilayer),'t23'))
           stheta = theta(ilayer);
       elseif(strcmp(stackings(ilayer),'r'))
           stheta = theta(ilayer);
       end
       str_stheta = num2str(stheta);
       filename = join([filename,str_stheta],"_theta=");
   end
   str_m = num2str(m);
   str_n = num2str(n);
   filename = join([filename,str_m],"_m=");
   filename = join([filename,str_n],"_n=");

   if (write_cart)
      if (strcmp(plot_f,'xyz'))
         write_xyz(nat,filename,alatL,id,pos)
         %write_full_xyz(nat,filename,alatL,id1,id2,id3,id4,pos1,pos2,pos3,pos4)
      elseif (strcmp(plot_f,'xsf'))
         write_xsf(nat,filename,alatL,z,id,num_at,atomlabels,pos,nlayers)
      end
   end
   
   if(write_positions)
       write_tb(nat,a,filename,alatL,alat1',theta,id,pos,nlayers);
   end
   
   % HEADER of geometry file for lammps
   if (write_lammps_input)
      write_lammps(nat,filename,alatL,z,id,3*nlayers,mass,atomlabels,pos,nlayers);
   end

end


function y1 = compdist(pos,pos2,nat,alatL,d1,d2,plot)
    n = 0;
    supercell = zeros(9*size(pos,1),2);
    supercell2 = zeros(9*size(pos,1),2);
    talat = alatL';
    for ni = -1 : 1
        for nj = -1:1
            for j = 1:nat/2
                n = n + 1;
                xshift = pos(j,1) + ni*talat(1,1) + nj*talat(1,2);
                yshift = pos(j,2) + ni*talat(2,1) + nj*talat(2,2);
                supercell(n,:) = [xshift,yshift];
                xshift2 = pos2(j,1) + ni*talat(1,1) + nj*talat(1,2);
                yshift2 = pos2(j,2) + ni*talat(2,1) + nj*talat(2,2);
                supercell2(n,:) = [xshift2,yshift2];
            end
        end
    end
      
    superd1 = zeros(9*size(nat/2,1),1);
    superd2 = zeros(9*size(nat/2,1),1);
    for i = 1 : 9
        superd1((i-1)*nat/2+1:i*nat/2) = d1(:);
        superd2((i-1)*nat/2+1:i*nat/2) = d2(:);
    end
    x = linspace(1.2*min(pos(:,1)),1.2*max(pos(:,1)),300);
    y = linspace(1.2*min(pos(:,2)),1.2*max(pos(:,2)),300);
    [X,Y] = meshgrid(x,y);
    f1 = scatteredInterpolant(supercell(:,1),supercell(:,2),superd1,'natural');
    Z1 = f1(X,Y);
    f2 = scatteredInterpolant(supercell2(:,1),supercell2(:,2),superd2,'natural');
    Z2 = f2(X,Y);
    if (plot)
        figure
        contourf(X,Y,Z2-Z1)
        figure
        surf(X,Y,Z2-Z1)
        shading interp
    end
    y1 = superd1;
end
