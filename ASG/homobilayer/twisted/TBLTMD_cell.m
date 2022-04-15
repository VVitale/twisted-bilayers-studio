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
function [nat,pos,pos2] = TBLTMD_cell(n, a, z, idMet, idChal, mass1, mass2, num_at1, num_at2, plot_f, kpath, nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase)
   % Setting parameters 
   m=n+1;
   
   theta = acos((n^2 + 4*n*m + m^2)/(2*(n^2 + n*m + m^2)));
   
   filename = "tBLTMD";
   filename = join([filename idMet],"_");
   filename = join([filename idChal "2"],'');

   % total number of atoms in Moire cell
   tot_nat= 6*(n^2 + n*m + m^2);
      
   % Compute Moire supercell  
   st=sin(theta);
   ct=cos(theta);
   u=0.6226;
   c=2*z;%12.2985;
   Rtheta = [ct -st; st ct];
   % Lattice vectors of top and bottom layer  
   ualat1 = a*[sqrt(3)/2, sqrt(3)/2;  -0.5,0.5];
   Rphi = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)];
   alat1 = Rphi*ualat1;
   alat2 = Rtheta*alat1;
   % Lattice vectors of Moire cell
   alatL = [n, m; -m, n+m]*alat1'
   ialatL = inv(alatL');
   
   % Generate reciprocal lattice vectors
   recL = 2*pi*inv(alatL');
   irecL = inv(recL);
   trecL = recL';
   ralat1 = 2*pi*inv(alat1');
   ralat2 = 2*pi*inv(alat2');

   % Print k-point grid/path to file
   theta = theta*180/pi;
   str_theta = num2str(theta);
   str_m = num2str(m);
   str_n = num2str(n);
   if (write_kpath)
      msg = ['Generating kpoint grid for moirè lattice with a twist angle of ',str_theta,' ...'];
      disp(msg)
      disp(' ')
      write_kpts(nk,ralat1,recL,str_theta,str_m,str_n,kpath,1);
      disp(' ')
      disp('... done')
   end
   
   % number of supercell for the search
   period = norm(alatL,1);
   tnat = floor(period/a);
   pos = zeros(tot_nat/2,3);
   pos2 = zeros(tot_nat/2,3);
   pos(:) = -99999;
   pos2(:) = -99999;
   id1(1,1:tot_nat/2)= "X";
   id2(1,1:tot_nat/2)= "X";
   delta =  1*10^(-5);
   
   disp(' ')
   msg = ['Generating Cartesian atomic coordinates for moirè lattice with a twist angle of ',str_theta,' ...'];
   disp(msg)

   k1=1;
   k2=1;
   for nj = -3*tnat:3*tnat
       for ni = -3*tnat:3*tnat
           r1 = ni*alat1(:,1) + nj*alat1(:,2);
           r2 = r1 + 1/3*alat1(:,1) + 1/3*alat1(:,2);
           r3 = r2;
           l1 = ialatL*r1;
           l2 = ialatL*r2;
           l3 = ialatL*r3;
           if (l1(1) < (1-delta) && l1(1) >= -delta && l1(2) < (1-delta) && l1(2) >= -delta) 
               pos(k1,1)=r1(1);
               pos(k1,2)=r1(2);
               pos(k1,3)=0.25*c;
               id1(k1) = idMet;
               k1 = k1 + 1;
           end
           if (l2(1) < (1-delta) && l2(1) >= -delta && l2(2) < (1-delta) && l2(2) >= -delta)
               pos(k1,1)=r2(1);
               pos(k1,2)=r2(2);
               pos(k1,3)= -(0.5-u)*c;
               id1(k1) = idChal;
               k1 = k1 + 1;
           end
           if (l3(1) < (1-delta) && l3(1) >= -delta && l3(2) < (1-delta) && l3(2) >= -delta) 
               pos(k1,1)=r3(1);
               pos(k1,2)=r3(2);
               pos(k1,3)=(1-u)*c;
               id1(k1) = idChal;
               k1 = k1 + 1;
           end
           if(phase=='AA')
              r1r = Rtheta*r1;
              r2r = Rtheta*r2;
              r3r = Rtheta*r3;
           elseif(phase=='AB')
              r1r = Rtheta*(-r1+1/3*alat1(:,1)+1/3*alat1(:,2));
              r2r = Rtheta*(-r2+1/3*alat1(:,1)+1/3*alat1(:,2));%[r2 + 2/3*alat1(:,1) + 2/3*alat1(:,2)];
              r3r = Rtheta*(-r2+1/3*alat1(:,1)+1/3*alat1(:,2));%[r2 + 2/3*alat1(:,1) + 2/3*alat1(:,2)];
           end
           l1 = ialatL*r1r;
           l2 = ialatL*r2r;
           l3 = ialatL*r3r;
           if (l1(1) < (1-delta) && l1(1) >= -delta && l1(2) < (1-delta) && l1(2) >= -delta) 
               pos2(k2,1) = r1r(1);
               pos2(k2,2) = r1r(2);
               pos2(k2,3) = 0.75*c;
               id2(k2) = idMet;
               k2 = k2 + 1;
           end
           if (l2(1) < (1-delta) && l2(1) >= -delta && l2(2) < (1-delta) && l2(2) >= -delta)
               pos2(k2,1) = r2r(1);
               pos2(k2,2) = r2r(2);
               pos2(k2,3) = u*c;
               id2(k2) = idChal;
               k2 = k2 + 1;
           end
           if (l3(1) < (1-delta) && l3(1) >= -delta && l3(2) < (1-delta) && l3(2) >= -delta)
               pos2(k2,1) = r3r(1);
               pos2(k2,2) = r3r(2);
               pos2(k2,3) = (1.5-u)*c;
               id2(k2) = idChal;
               k2 = k2 + 1;
           end
       end
   end
   
   nat = k1 + k2 - 2;
   
   %Check number of atoms is 4*(n^2 + n*m + m^2)
   if (nat ~= 6*(n^2 + n*m + m^2))
       nat
       6*(n^2 + n*m + m^2)
       error('Wrong number of atoms');
   end
   disp('... done')
   
   %figure
   rlri1 = zeros(nat/2,1);
   rlri2 = zeros(nat/2,1);
   %[rlri1,rlri2] = compLRIgraphene(pos,pos2,nat,alatL,a,false,3);
   
   if(phase=='AA')
       stheta = theta;
   elseif(phase=='AB')
       stheta = 60 - theta;
   end
   str_stheta = num2str(stheta);
   filename = join([filename,str_stheta],"_theta=");
   filename = join([filename,str_m],"_m=");
   filename = join([filename,str_n],"_n=");

   dmin=4.34;
   dmax=4.58;
   delta = 0.5*(dmax-dmin);

   if (write_cart)
      if (plot_f=='xyz')
         write_xyz(nat,filename,alatL,id1,id2,pos,pos2)
         write_full_xyz(nat,filename,alatL,id1,id2,pos,pos2)         
      elseif (plot_f=='xsf')
         write_xsf(nat,filename,alatL,z,id1,id2,num_at1,num_at2,idMet,idChal,pos,pos2)
      end
   end
   
   % HEADER of xyz coordinates file for TB
   if (write_positions)
      write_tb(nat,filename,alatL,alat1',stheta,id1,id2,pos,pos2)
   end
   
   % HEADER of geometry file for lammps
   if (write_lammps_input)
      m(1) = mass1;
      m(2) = mass2;
      m(3) = mass1; 
      m(4) = mass2;
      write_lammps(nat,filename,alatL,z,id1,id2,4,m,idMet,idChal,pos,pos2);
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
