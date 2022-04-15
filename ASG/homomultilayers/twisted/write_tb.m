function write_tb(nat,a,filename,ualatL,ucell,theta,id,upos,nlayers)
   % HEADER of xyz coordinates file for TB
   tbfname = 'positions';
   tbfname = join([tbfname,filename,'dat'],".");
   fileID = fopen(tbfname,'w');
   disp(' ')
   msg=['Writing coordinates and moire lattice vectors into TB input file:',tbfname];
   disp(msg)
   alatL = ualatL;
   pos = upos;
   cell = ucell;
   
   nat_per_layer_vec = nat/nlayers*ones(1,nlayers);
   format1 = '';
   for ilayer = 1 : nlayers
   format1 = join([format1,' %i ']);
   end
   format1 = join([format1,'\n']);
   
   format2 = '';
   for ilayer = 1 : nlayers
   format2 = join([format2,' %2.4f ']);
   end
   format2 = join([format2,'\n']);
   
   fprintf(fileID,format1,nat_per_layer_vec);
   fprintf(fileID,format2,a*ones(1,nlayers));
   fprintf(fileID,format2,theta);
   fprintf(fileID,format2,ones(1,nlayers));
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(1,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(2,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(1,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(2,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
   
   for ilayer = 0 : nlayers - 1
       for ik = 1:nat/nlayers
           k1 = ik + ilayer*nat/nlayers;
           fprintf(fileID,'%s %i %4.6f %4.6f %4.6f\n',id(k1),ilayer,pos(k1,1),pos(k1,2),pos(k1,3));
       end
   end
   
   fclose(fileID);
end
