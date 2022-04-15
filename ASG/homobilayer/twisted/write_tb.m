function write_tb(nat,filename,ualatL,ucell,theta,id1,id2,upos,upos2)
   % HEADER of xyz coordinates file for TB
   tbfname = 'positions';
   tbfname = join([tbfname,filename,'dat'],".");
   fileID = fopen(tbfname,'w');
   disp(' ')
   msg=['Writing coordinates and moire lattice vectors into TB input file:',tbfname];
   disp(msg)
   alatL = ualatL;
   pos = upos;
   pos2 = upos2;
   cell = ucell;
   
   fprintf(fileID,'%i %i\n',nat/2,nat/2);
   fprintf(fileID,'%2.4f\n',theta);
   fprintf(fileID,'%2.2f\n',1.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(1,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',cell(2,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(1,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',alatL(2,1:2), 0.0);
   fprintf(fileID,'%8.8f   %8.8f   %8.8f\n',0.0, 0.0,40);
   
   for k1 = 1:nat/2
      fprintf(fileID,'%s %i %4.6f %4.6f %4.6f\n',id1(k1),1,pos(k1,1),pos(k1,2),pos(k1,3));
   end
   
   for k1 = 1:nat/2
      fprintf(fileID,'%s %i %4.6f %4.6f %4.6f\n',id2(k1),2,pos2(k1,1),pos2(k1,2),pos2(k1,3));
   end
   
   fclose(fileID);
end
