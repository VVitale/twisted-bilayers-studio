function write_xyz(nat,filename,ualatL,id1,id2,upos,upos2)
   xyzfname = join([filename,"xyz"],".");
   fileID = fopen(xyzfname,'w');
   
   disp(' ')
   msg=['Writing xyz coordinates into file:',xyzfname];
   disp(msg)
   phi=acos(ualatL(1,1)/norm(ualatL(1,:)));
   cf=cos(-phi);
   sf=sin(-phi);
   Rphi=[cf -sf 0; sf cf 0; 0 0 1];
   alatL = ualatL*Rphi(1:2,1:2)';
   pos = upos*Rphi';
   pos2 = upos2*Rphi';
   % HEADER of XYZ file
   fprintf(fileID,'%i\n',nat);
   fprintf(fileID,'%s\n','');
   % Bottom layer first
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),pos(k1,1),pos(k1,2),pos(k1,3));
   end
   % Then top layer
   for k2 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
   end
   fclose(fileID);
end
