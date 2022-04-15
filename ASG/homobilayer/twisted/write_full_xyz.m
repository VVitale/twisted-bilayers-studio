function write_xyz(nat,filename,ualat,id1,id2,pos,pos2)
   xyzfname = join([filename,"_full",".xyz"]);
   fileID = fopen(xyzfname,'w');
   
   disp(' ')
   msg=['Writing xyz coordinates into file:',xyzfname];
   disp(msg)
   % HEADER of XYZ file
   fprintf(fileID,'%i\n',4*nat);
   fprintf(fileID,'%s\n','');
   % Bottom layer first
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),pos(k1,1),pos(k1,2),pos(k1,3));
   end
   % Then top layer
   for k2 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
   end
   shift = [ualat(1,:),0.0];
   tpos = pos + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos(k1,1),tpos(k1,2),tpos(k1,3));
   end
   tpos2 = pos2 + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos2(k1,1),tpos2(k1,2),tpos2(k1,3));
   end
   shift = [ualat(2,:),0.0];
   tpos = pos + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos(k1,1),tpos(k1,2),tpos(k1,3));
   end
   tpos2 = pos2 + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos2(k1,1),tpos2(k1,2),tpos2(k1,3));
   end
   shift = [ualat(1,:) + ualat(2,:),0.0];
   tpos = pos + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos(k1,1),tpos(k1,2),tpos(k1,3));
   end
   tpos2 = pos2 + shift.*ones(nat/2,1);
   for k1 = 1:nat/2
       fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),tpos2(k1,1),tpos2(k1,2),tpos2(k1,3));
   end
  % pos = pos*Rphi';
  % for k1 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),pos(k1,1),pos(k1,2),pos(k1,3));
  % end
  % pos = pos*Rphi';
  % for k1 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id1(k1),pos(k1,1),pos(k1,2),pos(k1,3));
  % end
  % pos2 = pos2*Rphi';
  % for k2 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
  % end
  %    pos2 = pos2*Rphi';
  % for k2 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
  % end
  %    pos2 = pos2*Rphi';
  % for k2 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
  % end
  %    pos2 = pos2*Rphi';
  % for k2 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
  % end
  %    pos2 = pos2*Rphi';
  % for k2 = 1:nat/2
  %     fprintf(fileID,'%s %4.6f %4.6f %4.6f\n',id2(k2),pos2(k2,1),pos2(k2,2),pos2(k2,3));
  % end
   fclose(fileID);
end
