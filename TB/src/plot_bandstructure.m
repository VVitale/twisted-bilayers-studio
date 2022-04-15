%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot bandstructure
% Adapted from:
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot bands

fprintf('--> Plotting bandstructure on Gamma-M-K-Gamma ... ')
line([scale_axis(knum),scale_axis(knum)],[-100,100],'Color','k');
hold on
line([scale_axis(2*knum),scale_axis(2*knum)],[-100,100],'Color','k');


if(reduced_workspace)
   plot(scale_axis,tb_bands,'ko','LineWidth',1.0,'MarkerSize',1.5);
else
   plot(scale_axis,tb_bands,'r','LineWidth',2.0);
end

%ylabel('Band Energy (eV)');
axis([-inf,inf,-6,5]);
box on;

ax = gca;
set(gca,'XTick',[0,scale_axis(knum),scale_axis(2*knum),1]);
set(gca,'YTick',[-6,-4,-2,0,2,4,6]);
if(~miniBZ)
    set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
else
    set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
end
set(gca,'YTickLabel',{'-6','-4','-2','0','2','4','6'});
set(gca,'FontSize',26);
set(gca,'Fontname','Times New Roman');
ylabel('Energy (eV)');
fprintf('done\n\n')
fprintf('--> Saving bandstructure in %s ... ',join([outfname,'.png']))
saveas(gca,join([outfname,'.png']))
fprintf('done\n\n')
%
if (twisted)
   figure
   line([scale_axis(knum),scale_axis(knum)],[-100,100],'Color','k');
   hold on
   line([scale_axis(2*knum),scale_axis(2*knum)],[-100,100],'Color','k');
   
   plot(scale_axis,tb_bands,'ko','LineWidth',1.0,'MarkerSize',2);
   box on;
   ax = gca;
   axis([-inf,inf,-0.85,-0.75]);
   set(gca,'XTick',[0,scale_axis(knum),scale_axis(2*knum),1]);
   set(gca,'YTick',[-0.85,-0.80,-0.75]);
   if(~miniBZ)
       set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
   else
       set(gca,'XTickLabel',{'\Gamma','K','M','\Gamma'});
   end    
   set(gca,'YTickLabel',{'-0.85','-0.80','-0.75'});
   set(gca,'FontSize',26);
   set(gca,'Fontname','Times New Roman');
   ylabel('Energy (eV)');
   
   saveas(gca,join([outfname,'_magnification.png']))
end
% 
% axis([-inf,inf,-0.1,0.1]);
% set(gca,'XTick',[0,scale_axis(knum),scale_axis(2*knum),1]);
% set(gca,'YTick',[-2,-1,0,1,2]);
% set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
% set(gca,'YTickLabel',{'-2','-1','0','1','2'});
% 
% saveas(gca,join([outfname,'_magnification_2.png']))
