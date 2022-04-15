function write_xsf(nat1,nat2,filename,alatL,c,at_num1,at_num2,pos,pos2,str_m,str_n,str_theta)

            xsffname = join([filename,str_theta],"_theta=");
            xsffname = join([xsffname,str_m],"_m=");
            xsffname = join([xsffname,str_n],"_n=");

            xsffname = join([xsffname "xsf"],".");
            fileID = fopen(xsffname,'w');
            disp(' ')
            msg=['Writing xsf coordinates into output file:',xsffname];
            disp(msg)
            % HEADER of XSF file
            fprintf(fileID,' %s\n','CRYSTAL');
            fprintf(fileID,' %s\n','PRIMVEC');
            fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',alatL(1,1), alatL(1,2), 0.0);
            fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',alatL(2,1), alatL(2,2), 0.0); 
            fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',0.0, 0.0, 5*c);
            fprintf(fileID,' %s\n','CONVVEC');
            fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',alatL(1,1), alatL(1,2), 0.0);
            fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',alatL(2,1), alatL(2,2), 0.0); 
            fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',0.0, 0.0, 5*c);
            fprintf(fileID,' %s\n','PRIMCOORD');
            fprintf(fileID,'\t%i\t%i\n',nat1+nat2,1);
            % Bottom layer first
            for k = 1:nat1
                    fprintf(fileID,' %i    %4.6f    %4.6f    %4.6f\n',at_num1(k),pos(k,1),pos(k,2),pos(k,3));
            end

            % Top layer
            for k = 1:nat2
                    fprintf(fileID,' %i   %4.6f    %4.6f    %4.6f\n',at_num2(k),pos2(k,1),pos2(k,2),pos2(k,3));
            end         
            fclose(fileID);
end
