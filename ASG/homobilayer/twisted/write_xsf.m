function write_xsf(nat,filename,ualatL,z,id1,id2,num_at1,num_at2,c_at1,c_at2,upos,upos2,d1,d2)
      xsffname = join([filename,"xsf"],".");
      fileID = fopen(xsffname,'w');
      disp(' ')
      msg=['Writing xsf coordinates into file:',xsffname];
      disp(msg)
      phi=acos(ualatL(1,1)/norm(ualatL(1,:)));
      cf=cos(-phi);
      sf=sin(-phi);
      Rphi=[cf -sf 0; sf cf 0; 0 0 1];
      alatL = ualatL*Rphi(1:2,1:2)';
      pos = upos*Rphi';
      pos2 = upos2*Rphi';
      % HEADER of XSF file
      fprintf(fileID,' %s\n','CRYSTAL');
      fprintf(fileID,' %s\n','PRIMVEC');
      fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',alatL(1,1), alatL(1,2), 0.0);
      fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',alatL(2,1), alatL(2,2), 0.0); 
      fprintf(fileID,' % 12.10f % 12.10f % 12.10f\n',0.0, 0.0, z^2 + 2);
      fprintf(fileID,' %s\n','CONVVEC');
      fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',alatL(1,1), alatL(1,2), 0.0);
      fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',alatL(2,1), alatL(2,2), 0.0); 
      fprintf(fileID,' %4.10f    %4.10f    %4.10f\n',0.0, 0.0, z^2 + 2);
      fprintf(fileID,' %s\n','PRIMCOORD');
      fprintf(fileID,'\t%i\t%i\n',nat,1);
      % Bottom layer first
      for k = 1:nat/2
          if (id1(k) == c_at1)
             idat = num_at1;
          else
             idat = num_at2;
          end
          fprintf(fileID,' %i    %4.6f    %4.6f    %4.6f\n',idat,pos(k,1),pos(k,2),pos(k,3));
      end

      % Top layer
      for k = 1:nat/2
          if (id2(k) == c_at1)
             idat = num_at1;
          else
             idat = num_at2;
          end
              fprintf(fileID,' %i   %4.6f    %4.6f    %4.6f\n',idat,pos2(k,1),pos2(k,2),pos2(k,3));
      end         
      fclose(fileID);

