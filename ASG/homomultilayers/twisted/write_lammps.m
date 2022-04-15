function write_lammps(nat,filename,alatL,z,id,ntype,m,atomlabels,pos,nlayers)
   % HEADER of geometry file for lammps
   lammpsfname = 'lammps_positions';
   lammpsfname = join([lammpsfname,filename,'dat'],".");
   additionalfname = join([filename,'rotation','dat'],".");
   fileID = fopen(lammpsfname,'w');
   fileID2 = fopen(additionalfname,'w');
   disp(' ')
   msg=['Writing coordinates, moire lattice info and masses into lammps input file:',lammpsfname];
   disp(msg)
   phi=acos(alatL(1,1)/norm(alatL(1,:)));
   cf=cos(-phi);
   sf=sin(-phi);
   Rphi=[cf -sf 0; sf cf 0; 0 0 1];
   rot_alatL = Rphi(1:2,1:2)*alatL';
   rot_pos = pos*Rphi';
   fprintf(fileID2,'%8.8f\n',phi*180/pi);
   fprintf(fileID2,'%8.8f %8.8f %8.8f\n',rot_alatL(:,1));
   fprintf(fileID2,'%8.8f %8.8f %8.8f\n',rot_alatL(:,2));
   fprintf(fileID,'%s\n','LAMMPS data file from TBLG_cell.m');
   fprintf(fileID,'\n');
   fprintf(fileID,'%i   %s\n',nat,'atoms');
   fprintf(fileID,'\n');
   fprintf(fileID,'%i   %s\n',ntype,'atom types');
   fprintf(fileID,'\n');
   fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0, rot_alatL(1,1), 'xlo  xhi');
   fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0, rot_alatL(2,2), 'ylo  yhi');
   fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0,100*z,'zlo  zhi');
   fprintf(fileID,'%8.8f   %8.8f   %8.8f   %s\n', rot_alatL(1,2), 0.0, 0.0, 'xy  xz  yz');
   fprintf(fileID,'\n');
   fprintf(fileID,'%s\n','Masses');
   fprintf(fileID,'\n');
   for i = 1 : ntype
      fprintf(fileID,'%i   %2.2f\n',i,m(i));
   end
   fprintf(fileID,'\n');
   fprintf(fileID,'%s\n','Atoms');
   fprintf(fileID,'\n');
   for ilayer = 0 : nlayers - 1
       for ik = 1:nat/nlayers
           k1 = ik+ilayer*nat/nlayers;
           if (strcmp(id(k1),atomlabels{3*ilayer+1}))
               fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',k1,1+3*ilayer,rot_pos(k1,1),rot_pos(k1,2),rot_pos(k1,3));
           elseif (strcmp(id(k1),atomlabels{3*ilayer+2}))
               fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',k1,2+3*ilayer,rot_pos(k1,1),rot_pos(k1,2),rot_pos(k1,3));
           elseif (strcmp(id(k1),atomlabels{3*ilayer+3}))
               fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',k1,3+3*ilayer,rot_pos(k1,1),rot_pos(k1,2),rot_pos(k1,3));
           end
       end
   end
   fclose(fileID);
   fclose(fileID2);
end
