function write_xyz(nat1,nat2,filename,id1,id2,pos,pos2,str_m,str_n,str_theta)
   filename = join([filename,str_theta],"_theta=");
   filename = join([filename,str_m],"_m=");
   filename = join([filename,str_n],"_n=");

   xyzfname = join([filename,"xyz"],".");
   fileID = fopen(xyzfname,'w');
   
   disp(' ')
   msg=['Writing xyz coordinates into file:',xyzfname];
   disp(msg)
   % HEADER of XYZ file
   fprintf(fileID,'%i\n',nat1+nat2);
   fprintf(fileID,'%s\n','');
   % Bottom layer first
   for k1 = 1:nat1
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),pos(k1,1),pos(k1,2),pos(k1,3));
   end
   % Then top layer
   for k2 = 1:nat2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
   end
   fclose(fileID);
end
