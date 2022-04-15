function [ax1,ax2,ax3,ax4] = plot_corrugation_quadlayer(class,system,a10,a20,a30)
p = a10/a20;
f = figure;
set(f,'Position',[500.5 500.5 800 600]);
t = tiledlayout(4,1,'TileSpacing','compact');
ax1=nexttile;
set(ax1,'FontSize',30);
set(ax1,'xticklabel',[]);
set(ax1,'LineWidth',4);
ax2=nexttile;
set(ax2,'FontSize',30);
set(ax2,'xticklabel',[]);
set(ax2,'LineWidth',4);
ax3=nexttile;
set(ax3,'FontSize',30);
set(ax3,'xticklabel',[]);
set(ax3,'LineWidth',4);
ax4=nexttile;
set(ax4,'FontSize',30);
set(ax4,'LineWidth',4);
legend('off');
hold(ax1,'on');hold(ax2,'on');hold(ax3,'on');hold(ax4,'on');
    c1 = [152, 78, 163]/255;
    c2 = [228, 26, 28]/255;
    c3 = [77, 175, 74]/255;
    c4 = [55, 126, 184]/255;
    c7 = [166, 206, 227]/255; %lightblue
    c5 = [178, 223, 138]/255; %lightgreen
    c6 = [253, 191, 111]/255; %lightorange
    c9 = [255, 127, 0]/255;   %orange
    c8 = [153, 153, 153]/255; %grey
    c10 = [247, 129, 191]/255; %pink
if(class=="homos")
%      [~,~,s3,t3] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BMX_t2H_BMX',c3);
%      [~,~,s4,t4] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BMX_t2H_BXM',c7);
%      [~,~,s5,t5] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BXM_t2H_BXM',c9);
%      [~,~,s1,t1] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_t2H_BMX',c1);
%      [~,~,s2,t2] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_t2H_BXM',c2);
%      [~,~,s2,t2] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_t2H_2H',c10);

     [~,~,s3,t3] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BMX_tAA_BMX',c1);
     [~,~,s4,t4] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BMX_tAA_BXM',c2);
     [~,~,s5,t5] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_BXM_tAA_BXM',c3);
     [~,~,s1,t1] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_tAA_BMX',c6);
     [~,~,s2,t2] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_tAA_BXM',c10);
     [~,~,s2,t2] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_2H_tAA_2H',c8);

%      [~,~,s3,t3] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_tAA_tAA_tAA',c1);
%      [~,~,s4,t4] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_t2H_tAA_tAA',c2);
%      [~,~,s5,t5] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,2,'homos','positions_t2H_tAA_t2H',c3);

    axis(ax1,[0 2 9.1 9.5]);
    axis(ax2,[0 2 2.9 3.3]);
    axis(ax3,[0 2 -3.3 -2.9]);
    axis(ax4,[0 2 -9.5 -9.1]);
    
%     legend(ax1,["B$^{M/X}$/t2H/B$^{M/X}$";...
%         "B$^{M/X}$/t2H/B$^{X/M}$";"B$^{X/M}$/t2H/B$^{X/M}$";...
%         "2H/t2H/B$^{M/X}$";"2H/t2H/B$^{X/M}$";"2H/t2H/2H"],'Location','northoutside',...
%         'FontSize',22,'Orientation','horizontal','NumColumns',3,...
%         'Interpreter','latex');

    legend(ax1,["B$^{M/X}$/tAA/B$^{M/X}$";...
        "B$^{M/X}$/tAA/B$^{X/M}$";"B$^{X/M}$/tAA/B$^{X/M}$";...
        "2H/tAA/B$^{M/X}$";"2H/tAA/B$^{X/M}$";"2H/tAA/2H"],'Location','northoutside',...
        'FontSize',22,'Orientation','horizontal','NumColumns',3,...
        'Interpreter','latex');

%     legend(ax1,["tAA/tAA/tAA";...
%         "t2H/tAA/tAA";"t2H/tAA/t2H"],'Location','northoutside',...
%         'FontSize',22,'Orientation','horizontal','NumColumns',3,...
%         'Interpreter','latex');


    %legend(ax1,['1x1';'2x2';'3x3'],'Location','northwest','FontSize',16);
    xlabel(ax4,'$|\mathbf{s}|/|\mathbf{t}_1+\mathbf{t}_2|$','FontSize',32,'Interpreter','latex');


elseif(class=="heteros")
%       [~,~,s1,t1]=analyse_structure_supercell(ax1,ax2,2,class,'positions.MoSe2_on_MoS2_theta=4.531_m=4_n=14_rotated',[255 193 37]/255);
% %      [~,~,s2,t2]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t4',[255 193 37]/255);
%       [~,~,s3,t3]=analyse_structure_supercell(ax1,ax2,2,class,'positions.MoSe2_on_MoS2_theta=5.397_m=4_n=12_rotated',[155 35 53]/255);
% %      [~,~,s4,t4]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t2',[255 193 37]/255);
%       [~,~,s5,t5]=analyse_structure_supercell(ax1,ax2,2,class,'positions.MoSe2_on_MoS2_theta=7.89_m=5_n=4_rotated',[67 110 238]/255);

      [~,~,s1,t1]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t1_new',[255 193 37]/255);
      [~,~,s2,t2]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t2_new',[155 35 53]/255);%[255 193 37]/255);
      [~,~,s3,t3]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t3_new',[67 110 238]/255);
      
      %[~,~,s4,t4]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t3',[32 220 170]/255);
      %[~,~,s5,t5]=analyse_structure_supercell(ax1,ax2,2,class,'positions_t5',[180 42 248]/255);
% [~,~,s1,t1]=analyse_structure_supercell(ax1,ax2,2,class,'pt1',[0.5 1 0.5]);
% [~,~,s2,t2]=analyse_structure_supercell(ax1,ax2,2,class,'pt2',[1 0 1]);
% [~,~,s3,t3]=analyse_structure_supercell(ax1,ax2,2,class,'pt3',[1 0 1]);
% [~,~,s4,t4]=analyse_structure_supercell(ax1,ax2,2,class,'pt4',[1 0 1]);
% [~,~,s5,t5]=analyse_structure_supercell(ax1,ax2,2,class,'pt5',[1 0 1]);
% [~,~,s6,t6]=analyse_structure_supercell(ax1,ax2,2,class,'pt6',[1 0 1]);
% [~,~,s7,t7]=analyse_structure_supercell(ax1,ax2,2,class,'pt7',[1 0 1]);
% [~,~,s8,t8]=analyse_structure_supercell(ax1,ax2,2,class,'pt8',[1 0 1]);
% [~,~,s9,t9]=analyse_structure_supercell(ax1,ax2,2,class,'pt9',[1 0 1]);
% [~,~,s10,t10]=analyse_structure_supercell(ax1,ax2,2,class,'pt10',[1 0 1]);
%[~,~,s10,t10]=analyse_structure(ax1,ax2,2,class,'pt10',[rand(1) rand(1) rand(1)]);
%[~,~,s11,t11]=analyse_structure(ax1,ax2,2,class,'pt11',[rand(1) rand(1) rand(1)]);

%    [~,~,s1,t1]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.030',[rand(1) rand(1) rand(1)]);
%    [~,~,s2,t2]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.025',[(2-rand(1))/3 (2-rand(1))/3 (2-rand(1))/3]);
%    [~,~,s,t]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.020',[1 0.75 0]);
%    [~,~,s4,t4]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.015',[1 0.5 0]);
%    [~,~,s5,t5]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.010',[1 0.25 0]);
%    [~,~,s6,t6]=analyse_structure_supercell(ax1,ax2,2,class,'positions_-.005',[1 0 0]);
%    [~,~,s7,t7]=analyse_structure_supercell(ax1,ax2,2,class,'positions_0',[0 0 0]);
%    [~,~,s8,t8]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.005',[0 1 1]);
%    [~,~,s9,t9]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.010',[0 0.75 1]);
%    [~,~,s10,t10]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.015',[0 0.5 1]);
%    [~,~,s11,t11]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.020',[0 0.25 1]);
%    [~,~,s12,t12]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.025',[(2+rand(1))/3 (2+rand(1))/3 (2+rand(1))/3]);
%    [~,~,s13,t13]=analyse_structure_supercell(ax1,ax2,2,class,'positions_.030',[rand(1) rand(1) rand(1)]);
    %analyse_structure(ax1,ax2,2,class,'positions_.35',[rand(1) rand(1) rand(1)])
    %analyse_structure(ax1,ax2,2,class,'positions_.40',[rand(1) rand(1) rand(1)])
    %analyse_structure(ax1,ax2,2,class,'positions_.45',[rand(1) rand(1) rand(1)])
    %analyse_structure(ax1,ax2,2,class,'positions_.50',[rand(1) rand(1) rand(1)])
    
    %axis(ax1,[-inf inf 0 3]);
    %axis([-inf inf -4 -1]);
    %lgd=legend(ax1,[join([num2str(2.6),'°, ',num2str(((s1/100+1)/p-1)*100,'%2.3f\n'),"%"]);join([num2str(3.4),'°, ',num2str(((s2/100+1)/p-1)*100,'%2.3f\n'),"%"]);join([num2str(7.9),'°, ',num2str(((s3/100+1)/p-1)*100,'%2.3f\n'),"%"])],'Location','northwest','FontSize',28);
    %lgd=legend(ax1,[join([num2str(t1,'%2.1f\n'),' , ',num2str(s1,'%2.1f\n')]);join([num2str(t3,'%2.1f\n'),' , ',num2str(s3,'%2.1f\n')]);join([num2str(t5,'%2.1f\n'),' , ',num2str(s5,'%2.1f\n')])],'Location','northwest','FontSize',28);
    %lgd=legend(ax1,[join([num2str(t1,'%2.1f\n'),'°(',num2str(s1,'%2.1f\n'),' %)']);join([num2str(t2,'%2.1f\n'),'°(',num2str(s2,'%2.1f\n'),' %)']);join([num2str(t3,'%2.1f\n'),'°(',num2str(s3,'%2.1f\n'),' %)']);join([num2str(t4,'%2.1f\n'),'°(',num2str(s4,'%2.1f\n'),' %)']);join([num2str(t5,'%2.1f\n'),'°(',num2str(s5,'%2.1f\n'),' %)'])],'Location','north','FontSize',22,'Orientation','horizontal','NumColumns',3);
    lgd=legend(ax1,[join([num2str(t1,'%2.1f\n'),'°(',num2str(s1,'%2.1f\n'),' %)']);join([num2str(t2,'%2.1f\n'),'°(',num2str(s2,'%2.1f\n'),' %)']);join([num2str(t3,'%2.1f\n'),'°(',num2str(s3,'%2.1f\n'),' %)'])],'Location','northoutside','FontSize',28,'Orientation','horizontal','NumColumns',3);
    %lgd.ItemTokenSize = [15,10];
    %lgd = legend(ax1,[join([num2str(-3.0,'%2.1f\n'),"%"," , ",num2str(s-3.0,'%2.1f\n')]);join([num2str(-2.5,'%2.1f\n'),"%"," , ",num2str(s-2.5,'%2.1f\n')]);join([num2str(-2.0,'%2.1f\n'),"%"," , ",num2str(s-2.0,'%2.1f\n'),"%"]); join([num2str(-1.5,'%2.1f\n'),"%"," , ",num2str(s-1.5,'%2.1f\n'),"%"]); join([num2str(-1.0,'%2.1f\n'),"%"," , ",num2str(s-1.0,'%2.1f\n'),"%"]); join([num2str(-0.5,'%2.1f\n'),"%"," , ",num2str(s-0.5,'%2.1f\n'),"%"]); join([num2str(0.0,'%2.1f\n'),"%"," , ",num2str(s,'%2.1f\n'),"%"]); join([num2str(+0.5,'%2.1f\n'),"%"," , ",num2str(s+0.5,'%2.1f\n'),"%"]); join([num2str(+1.0,'%2.1f\n'),"%"," , ",num2str(s+1,'%2.1f\n'),"%"]); join([num2str(+1.5,'%2.1f\n'),"%"," , ",num2str(s+1.5,'%2.1f\n'),"%"]); join([num2str(+2.0,'%2.1f\n'),"%"," , ",num2str(s+2,'%2.1f\n'),"%"]);join([num2str(2.5,'%2.1f\n'),"%"," , ",num2str(s+2.5,'%2.1f\n')]);join([num2str(3.0,'%2.1f\n'),"%"," , ",num2str(s+3.0,'%2.1f\n')])],'Location','northeastoutside');
    %lgd = legend(ax1,[join([num2str(-2.0,'%2.1f\n'),"%"," , ",num2str(s-2.0,'%2.1f\n'),"%"]); join([num2str(-1.5,'%2.1f\n'),"%"," , ",num2str(s-1.5,'%2.1f\n'),"%"]); join([num2str(-1.0,'%2.1f\n'),"%"," , ",num2str(s-1.0,'%2.1f\n'),"%"]); join([num2str(-0.5,'%2.1f\n'),"%"," , ",num2str(s-0.5,'%2.1f\n'),"%"]); join([num2str(0.0,'%2.1f\n'),"%"," , ",num2str(s,'%2.1f\n'),"%"]); join([num2str(+0.5,'%2.1f\n'),"%"," , ",num2str(s+0.5,'%2.1f\n'),"%"]); join([num2str(+1.0,'%2.1f\n'),"%"," , ",num2str(s+1,'%2.1f\n'),"%"]); join([num2str(+1.5,'%2.1f\n'),"%"," , ",num2str(s+1.5,'%2.1f\n'),"%"]); join([num2str(+2.0,'%2.1f\n'),"%"," , ",num2str(s+2,'%2.1f\n'),"%"])],'Location','northeastoutside');
    %lgd = legend(ax1,[join([num2str(-2.0,'%2.1f\n'),"%"," , ",num2str(s-2.0,'%2.1f\n'),"%"])],'Location','northeastoutside');
    %lgd=legend(ax1,["1x1 from 2x2 no\_relax";"1x1 from 2x2 relaxed"],'Location','northwest','FontSize',28);
   %lgd_txt=[join([num2str(t1,'%2.3f\n'),'°, ',num2str(((s1/100+1)/p-1)*100,'%2.3f\n'),"%"])];
        %join([num2str(t2,'%2.3f\n'),'°, ',num2str(((s2/100+1)/p-1)*100,'%2.3f\n'),"%"])
        %join([num2str(t3,'%2.3f\n'),'°, ',num2str(((s3/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t4,'%2.3f\n'),'°, ',num2str(((s4/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t5,'%2.3f\n'),'°, ',num2str(((s5/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t6,'%2.3f\n'),'°, ',num2str(((s6/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t7,'%2.3f\n'),'°, ',num2str(((s7/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t8,'%2.3f\n'),'°, ',num2str(((s8/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t9,'%2.3f\n'),'°, ',num2str(((s9/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t10,'%2.3f\n'),'°, ',num2str(((s10/100+1)/p-1)*100,'%2.3f\n'),"%"])
%         join([num2str(t11,'%2.3f\n'),'°, ',num2str(((s11/100+1)/p-1)*100,'%2.3f\n'),"%"])];
        %join([num2str(t12,'%2.3f\n'),'°, ',num2str(((s12/100+1)/p-1)*100,'%2.3f\n'),"%"])
        %join([num2str(t13,'%2.3f\n'),'°, ',num2str(((s13/100+1)/p-1)*100,'%2.3f\n'),"%"])];
    xlabel(ax2,"$|\mathbf{s}|/|\mathbf{t}_2'-\mathbf{t}_1'|$",'FontSize',32,'Interpreter','latex');

    %lgd=legend(ax1,lgd_txt,'Location','northeastoutside','FontSize',20);
    title(lgd,'   \theta(\epsilon)','FontSize',28,'FontWeight','bold')
    axis(ax1,[0 2 1.4 2])
    axis(ax2,[0 2 -2 -1.4])
    %axis(ax1,[0 1 3.1 3.7])
    %axis(ax2,[0 1 -3.7 -3.1])
    set(ax1,'Ytick',(0:0.5:5));
    set(ax2,'Ytick',(-5:0.5:0));
elseif(class=='all')


    [~,~,s1,t1]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_mos2_new',c1);
    [~,~,s2,t2]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_mose2_new',c2);
    [~,~,s3,t3]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_ws2_new',c3);
    [~,~,s4,t4]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_wse2_new',c4);
    %[~,~,s5,t5]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_mos2ws2_new',c5,3,2);
    %[~,~,s6,t6]=analyse_structure_supercell(ax1,ax2,2,'homos','positions_mose2wse2_new',c6,3,2);
    %[~,~,s7,t7]=analyse_structure_supercell(ax1,ax2,2,'heteros','positions_mos2mose2_new',c7,3,3);
    %[~,~,s8,t8]=analyse_structure_supercell(ax1,ax2,2,'heteros','positions_ws2mose2_new',c8,3,3);
    %[~,~,s9,t9]=analyse_structure_supercell(ax1,ax2,2,'heteros','positions_mos2wse2_new',c9,3,3);
    %[~,~,s10,t10]=analyse_structure_supercell(ax1,ax2,2,'heteros','positions_ws2wse2_new',c10,3,3);
    
    %xlabel(ax2,"$|\mathbf{s}|/|\mathbf{t}_2'-\mathbf{t}_1'|$",'FontSize',32,'Interpreter','latex');
    xlabel(ax2,'$|\mathbf{s}|/|\mathbf{t}_1+\mathbf{t}_2|$','FontSize',32,'Interpreter','latex');
    legend(ax1,'MoS$_2$/MoS$_2$','MoSe$_2$/MoSe$_2$','WS$_2$/WS$_2$','WSe$_2$/WSe$_2$','Location','Northoutside','FontSize',29,'Interpreter','Latex');
    %legend(ax1,'WS$_2$/MoS$_2$','WSe$_2$/MoSe$_2$','MoSe$_2$/MoS$_2$','MoSe$_2$/WS$_2$','WSe$_2$/MoS$_2$','WSe$_2$/WS$_2$','Location','Northoutside','FontSize',28,'Interpreter','Latex');
    %legend(ax1,'MoSe$_2$/MoS$_2$','MoSe$_2$/WS$_2$','WSe$_2$/MoS$_2$','WSe$_2$/WS$_2$','Location','NorthEastoutside','FontSize',28,'Interpreter','Latex');

    %lgd=legend(ax1,lgd_txt,'Location','northeastoutside','FontSize',20);
    %title(lgd,'   \theta(°) \sigma(%)','FontSize',28,'FontWeight','bold')
    %axis(ax1,[0 2 2.9 3.7])
    %axis(ax2,[0 2 -3.7 -2.9])
    axis(ax1,[0 2 1 1.8])
    axis(ax2,[0 2 -2.6 -1])
    set(ax1,'Ytick',(0:0.5:5));
    set(ax2,'Ytick',(-5:0.5:0));
end
set(ax1,'Box','on','Linewidth',4);
set(ax2,'Box','on','Linewidth',4);
set(ax3,'Box','on','Linewidth',4);
set(ax4,'Box','on','Linewidth',4);
a20/a10-1
t.Title.String = join([system]);
set(t.Title,'Interpreter','latex')
t.Title.FontSize = 40;
ylabel(t,'$\Delta$z(\AA)','FontSize',36,'Interpreter','latex');


end
