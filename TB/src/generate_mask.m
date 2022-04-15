function M = generate_mask(orbitals,tot_norbs,nlayer,flipped)
M = eye(tot_norbs);
I = ones(tot_norbs,1);
for ilayer = 1 : nlayer
    if(flipped(ilayer)==1)
       orby = findobj(orbitals,'Layer',ilayer,{'Rel_index',2,'-or','Rel_index',5,'-or','Rel_index',7,'-or','Rel_index',11});
       lorby = [orby(:).Ham_index];
       I(lorby) = -1;
    end
end
M = M.*I;
end  
