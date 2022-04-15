function write_pot(nat1,nat2,filename,id1,id2,id11,id12,id13,id21,id22,id23,pot_11,pot_12,pot_13,pot_21,pot_22,pot_23,pos,pos2,str_m,str_n,str_theta)
        filename = join([filename,str_theta],"_theta=");
        filename = join([filename,str_m],"_m=");
        filename = join([filename,str_n],"_n=");

        disp(' ')
        msg = ['Generating potential.dat'];
        disp(msg)
        disp(' ')

        potfname = 'potential';
        potfname = join([potfname,filename,'dat'],".");
        fileID = fopen(potfname,'w');
        % Bottom layer
        for k = 1:nat1
               if (id1(k) == id11)
                  fprintf(fileID,'%2.4f\n',pot_11);
               elseif (id1(k) == id12)
                  fprintf(fileID,'%2.4f\n',pot_12);
               elseif (id1(k) == id13)
                  fprintf(fileID,'%2.4f\n',pot_13);
               end
        end

        % Top layer
        for k = 1:nat2
               if (id2(k) == id21)
                  fprintf(fileID,'%2.4f\n',pot_21);
               elseif (id2(k) == id22)
                  fprintf(fileID,'%2.4f\n',pot_22);
               elseif (id2(k) == id23)
                  fprintf(fileID,'%2.4f\n',pot_23);
               end
        end
        fclose(fileID);
end
