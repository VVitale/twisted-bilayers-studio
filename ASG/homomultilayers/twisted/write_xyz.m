function write_xyz(nat,filename,ualatL,id,upos)
   xyzfname = join([filename,"xyz"],".");
   fileID = fopen(xyzfname,'w');
   
   disp(' ')
   msg=['Writing xyz coordinates into file:',xyzfname];
   disp(msg)
   phi=acos(ualatL(1,1)/norm(ualatL(1,:)));
   cf=cos(-phi);
   sf=sin(-phi);
   Rphi=[cf -sf 0; sf cf 0; 0 0 1];
   pos = upos*Rphi';
   % HEADER of XYZ file
   fprintf(fileID,'%i\n',nat);
   fprintf(fileID,'%s\n','');
   for k1 = 1:nat
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id(k1),pos(k1,1),pos(k1,2),pos(k1,3));
   end
   fclose(fileID);
end
