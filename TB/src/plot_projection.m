function y = plot_projection(orbitals,tot_norbs,knum_tot,knum,tb_bands,...
    tb_vecs,scale_axis,property,value,fname,stitle,end_color)

line([scale_axis(knum),scale_axis(knum)],[-100,100],'Color','k');
hold on
line([scale_axis(2*knum),scale_axis(2*knum)],[-100,100],'Color','k');

[~,~] = projection(orbitals,scale_axis,tb_bands,tb_vecs,fname,{property,value});

data = importdata(fname);
data_r = reshape(data(:,3),[knum_tot, tot_norbs]);
size(data_r)

for iband = 1 : tot_norbs
    h = surface([scale_axis';scale_axis'],...
        [tb_bands(iband,:);tb_bands(iband,:)],...
        [zeros(size(scale_axis'));zeros(size(scale_axis'))],...
        [data_r(:,iband)';data_r(:,iband)'],'FaceColor','none',...
        'EdgeColor','interp');
    set(h,'LineWidth',3)
    hold on
end
map = colorGradient([200 200 200]/255,end_color,100);
colormap(map)
colorbar
y = 0;
axis([-inf,inf,-6,5]);
box on;
ax = gca;
set(ax,'XTick',[0,scale_axis(knum),scale_axis(2*knum),1]);
set(ax,'YTick',[-6,-4,-2,0,2,4,6]);
set(ax,'XTickLabel',{'\Gamma','M','K','\Gamma'});
set(ax,'YTickLabel',{'-6','-4','-2','0','2','4','6'});
set(ax,'FontSize',26);
set(ax,'Fontname','Times New Roman');
ylabel('Energy (eV)');
title(stitle,'Interpreter','latex')
saveas(ax,join([fname,'.fig']))
end