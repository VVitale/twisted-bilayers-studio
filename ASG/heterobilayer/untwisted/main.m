% INPUT FILES
% Here we probably want a table of possible approximants
theta = 0.0;
n=25 ;
m=24 ;
a = 3.183;%/sqrt(3);
u = 0.621;
c = 12.32;
layer1='TMD'
layer2='TMD';

% Info for writing Cartesian coordinates and potential to file
write_positions = true;
write_potential = true;
write_cart      = true;
% Whether to XYZ or XSF format
plot_f = 'xyz';
write_lammps_input = false;

% Atomic symbol of each atomic species
id11 ="Mo";%'Mo';
id12 ="S";%'S';
id13 ="S";%'S';
id21 ="Mo";%'W';
id22 ="Se";%'Se';
id23 ="Se";%'Se';
% Atomic number of each species (used in XSF file)
at_num11 = 42 ;%42;
at_num12 = 16 ;%16;
at_num13 = 16 ;%16;
at_num21 = 42 ;%74;
at_num22 = 34 ;%34;
at_num23 = 34 ;%34;

% Mass of each atomic species in a.u. (used in LAMMPS input file)
mass_1 =95.94;
mass_2 =32.06;
mass_3 =95.94;%183.84;
mass_4 =78.96;

% Potential for each atomic species
pot_11 = 3.34;
pot_12 = -1.40;
pot_13 = -1.40;
pot_21 = 0.0;
pot_22 = 0.0;
pot_23 = 0.0;

write_kpath = true;
path = 'KGMKp';
init_kpts = 15;

% SOME PARAMETERS

str_theta = num2str(theta);
if(layer1=="TMD")
   f1 = join([id11,id12,"2"],"");
else
   f1 = join([id11,id12],"");
end
if(layer2=="TMD")
   f2 = join([id21,id22,"2"],"");
else
   f2 = join([id21,id22],"");
end
filename = join([f1,f2],"_on_");
filename = join([filename,str_theta],"_theta=");

if(layer1=="TMD" && layer2=="TMD")
   comb = 1;
elseif(layer1=="graphene" && layer2=="graphene")
   comb = 2;
elseif((layer1=="graphene" && layer2=="TMD"))
   comb = 3;
elseif (layer2=="graphene" && layer1=="TMD")
   comb = 4;
end

theta= theta*pi/180;
str_m = num2str(m);
str_n = num2str(n);
st=sin(theta);
ct=cos(theta);

Rtheta = [ct -st; st ct];
M = eye(2);
chi = m/n;
expansion = chi; %56/55;%ahBN/a;
M = M*(expansion);
alat1 = [1 0.5; 0 sqrt(3)/2]*a
alat2 = M*Rtheta*alat1

alatL = (eye(2) - inv(Rtheta)*inv(M))\alat1';
ialatL = inv(alatL');
intL1 = alatL/(alat1')
intL2 = alatL/(alat2')
if (any(abs(intL2 - round(intL2))>0.001))
    error('Wrong moirè lattice vectors')
end

tnat=100;
corr_nat = num_atoms(int64(m),int64(n),int64(comb));
pos=zeros(corr_nat,3);
pos2=zeros(corr_nat,3);
pos(:,:) = -99999;
pos2(:,:) = -99999;
id1(1,1:corr_nat) = "X";
id2(1,1:corr_nat) = "X";
% Array containing the atoms' numbers
at_num1(1,1:corr_nat) = 0;
at_num2(1,1:corr_nat) = 0;

delta = 1*10^(-5);

[nat,nat1,nat2,pos,pos2,id1,id2,at_num1,at_num2]=generate_structure(a,u,c,M,Rtheta,tnat,alat1,ialatL,delta,pos,pos2,id1,id2,id11,id12,id13,id21,id22,id23,at_num1,at_num2,at_num11,at_num12,at_num13,at_num21,at_num22,at_num23,layer1,layer2);

disp('done')
% %Check number of atoms is correct
% %if (nat ~= corr_nat)
% %    msg = ['Number of atoms is incorrect. Expected ',num2str(corr_nat),' found ',num2str(nat)];
% %    disp(msg)
% %    error('Error')
% %else
% %    msg = ['Number of atoms is correct : ', num2str(nat)];
% %    disp(msg)
% %end

% Generate reciprocal lattice
if (write_kpath)
   msg = ['Generating kpoint grid for moirè lattice ',num2str(theta*180/pi),' ...'];
   disp(msg)
   % Generate reciprocal lattice
   recL = 2*pi*inv(alatL');
   irecL = inv(recL);
   trecL = recL'; 
   ralat1 = 2*pi*inv(alat1');
   ralat2 = 2*pi*inv(alat2');
   write_kpts(init_kpts,ralat2,recL,str_theta,str_m,str_n,path,2)
   disp('done')
end

pos = pos(1:nat1,:);
pos2 = pos2(1:nat2,:);
id1 = id1(1,1:nat1);
id2 = id2(1,1:nat2);

%%%% WRITE POSITION FILE FOR TB_FREE.X %%%%
if (write_positions)
   write_pos(nat1,nat2,filename,alatL,c,id1,id2,pos,pos2,str_m,str_n,str_theta)
end

%%%% WRITE POTENTIAL FILE FOR TB_FREE.X %%%%
if (write_potential)
   write_pot(nat1,nat2,filename,id1,id2,id11,id12,id13,id21,id22,id23,pot_11,pot_12,pot_13,pot_21,pot_22,pot_23,pos,pos2,str_m,str_n,str_theta)
end
     
%%%% WRITE COORDINATES FILE %%%%
if(write_cart)
   if (plot_f=='xyz')
      write_xyz(nat1,nat2,filename,id1,id2,pos,pos2,str_m,str_n,str_theta)
   elseif (plot_f=='xsf')
      write_xsf(nat1,nat2,filename,alatL,c,at_num1,at_num2,pos,pos2,str_m,str_n,str_theta)
   end
end

%%%% WRITE LAMMPS INPUT FILE %%%%
if (write_lammps_input)
   m(1) = mass_1;
   m(2) = mass_2;
   m(3) = mass_3;
   m(4) = mass_4;

   write_lammps(nat1,nat2,filename,alatL,c,ntype,m,id1,id2,id11,id12,id13,id21,id22,id23,pos,pos2,str_m,str_n,str_theta)
end

function y = num_atoms(m,n,c)
    if ~isinteger(m)
        error('First imput must be an integer (m)')
    end
    if ~isinteger(n)
        error('Second imput must be an integer (n)')
    end
    if ~isinteger(c)
        error('Third imput must be an integer (c)')
    elseif c< -1 || c > 4 || c==0
        error('c must be equal to: 1 or 2 or 3 or 4')
    end
    
    if (c==1) 
        y = 3*(n^2) +3*(m^2);
    elseif (c==2)
        y = 2*(n^2) + 2*(m^2);
    elseif (c==3)
        y = 2*(n^2) + 3*(m^2);
    elseif (c==3)
        y = 3*(n^2) + 2*(m^2);
    end
    
end
