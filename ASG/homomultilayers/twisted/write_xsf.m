function write_xsf(nat,filename,ualatL,z,id,num_at,atomlabels,upos,nlayers)
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
      for ilayer = 0 : nlayers - 1
          for ik = 1:nat/nlayers
              k = ik+ilayer*nat/nlayers;
              if (strcmp(id(k),atomlabels{ilayer+1}))
                  idat = num_at(ilayer+1);
              elseif(strcmp(id(k),atomlabels{ilayer+2}))
                  idat = num_at(ilayer+2);
              elseif(strcmp(id(k),atomlabels{ilayer+3}))
                  idat = num_at(ilayer+3);
              end
              fprintf(fileID,' %i    %4.6f    %4.6f    %4.6f\n',idat,pos(k,1),pos(k,2),pos(k,3));
          end
      end
      fclose(fileID);
end

