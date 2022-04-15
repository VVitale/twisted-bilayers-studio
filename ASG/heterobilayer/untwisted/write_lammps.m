function write_lammps(nat1,nat2,filename,alatL,c,ntype,m,id1,id2,id11,id12,id13,id21,id22,id23,pos,pos2,str_m,str_n,str_theta)
        filename = join([filename,str_theta],"_theta=");
        filename = join([filename,str_m],"_m=");
        filename = join([filename,str_n],"_n=");
        lammpsfname = 'lammps_positions';
        lammpsfname = join([lammpsfname,filename,'dat'],".");
        fileID = fopen(lammpsfname,'w');
        % HEADER of positions.dat
        phi=acos(alatL(1,1)/norm(alatL(1,:)));
        cf=cos(-phi);
        sf=sin(-phi);
        Rphi=[cf -sf; sf cf];
        rot_alatL = Rphi*alatL';
        rot_pos = Rphi*pos(:,1:2)';
        rot_pos2 = Rphi*pos2(:,1:2)';
        fprintf(fileID,'%s\n','LAMMPS data file from TBLG_cell.m');
        fprintf(fileID,'\n');
        fprintf(fileID,'%i   %s\n',nat1+nat2,'atoms');
        fprintf(fileID,'\n');
        fprintf(fileID,'%i   %s\n',ntype,'atom types');
        fprintf(fileID,'\n');
        fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0, rot_alatL(1,1), 'xlo  xhi');
        fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0, rot_alatL(2,2), 'ylo  yhi');
        fprintf(fileID,'%8.8f   %8.8f   %s\n',0.0,10*c,'zlo  zhi');
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
        i = 1;
        for k = 1:nat1
                if (id1(k) == id11)
                   fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',k,1,rot_pos(1,k),rot_pos(2,k),pos(k,3));
                elseif (id1(k) == id12 || id1(k) == id13) 
                   fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',k,2,rot_pos(1,k),rot_pos(2,k),pos(k,3));
                end
        end
        
        i = 1;
        for k = 1:nat2
                if (id2(k) == id21)
                   fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',nat1+k,3,rot_pos2(1,k),rot_pos2(2,k),pos2(k,3));
                elseif (id2(k) == id22 || id2(k) == id23) 
                   fprintf(fileID,'%i %i %4.6f %4.6f %4.6f\n',nat1+k,4,rot_pos2(1,k),rot_pos2(2,k),pos2(k,3));
                end
        end
        
        fclose(fileID);
end
