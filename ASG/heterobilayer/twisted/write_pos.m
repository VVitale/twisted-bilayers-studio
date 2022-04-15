function write_pos(nat1,nat2,filename,alatL,theta,strain,cell,id1,id2,pos,pos2,str_m,str_n,str_theta)
        posfname = join([filename,str_theta],"_theta=");
        posfname = join([posfname,str_m],"_m=");
        posfname = join([posfname,str_n],"_n=");
        posfname = join(['positions',posfname,'dat'],".");
        fileID = fopen(posfname,'w');
        % HEADER of positions.dat
        fprintf(fileID,'%i %i\n',nat1,nat2);
        fprintf(fileID,'%2.4f\n',theta);
        fprintf(fileID,'%2.4f\n',strain);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(1,1:2), 0.0);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(2,1:2), 0.0);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(1,1:2), 0.0);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(2,1:2), 0.0);
        fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
        for k = 1:nat1
                %d1(k) = 1.0-delta*rlri1(k);
                fprintf(fileID,'%s %i %4.6f %4.6f %4.6f\n',id1(k),1,pos(k,1),pos(k,2),pos(k,3));
        end

        % Top layer
        for k = 1:nat2
                %d2(k) = dmin+delta*rlri2(k);
                fprintf(fileID,'%s %i %4.6f %4.6f %4.6f\n',id2(k),2,pos2(k,1),pos2(k,2),pos2(k,3));
        end
        fclose(fileID);
end
