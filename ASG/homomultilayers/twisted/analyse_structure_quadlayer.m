function [out1,out2,strain,theta] = analyse_structure_quadlayer(ax1,ax2,ax3,ax4,mode,class,inputfile,rgb)
data = importdata(join([inputfile,'.dat']),' ',10);
natoms = str2num(data.textdata{1});
lattice = str2num(data.textdata{2});
strain=(str2num(data.textdata{4})-1)*100
thetas=str2num(data.textdata{3});
theta = thetas(2);
nat1 = natoms(1);
nat2 = natoms(2);
nat3 = natoms(3);
nat4 = natoms(4);
pos1 = data.data(1:nat1,2:4);
pos2 = data.data(nat1+1:nat1+nat2,2:4);
pos3 = data.data(nat1+nat2+1:nat1+nat2+nat3,2:4);
pos4 = data.data(nat1+nat2+nat3+1:nat1+nat2+nat3+nat4,2:4);

cpos1 = mean(pos1(:,3)); %pos1(index1,3)
cpos2 = mean(pos2(:,3)); %pos2(index2,3)
cpos3 = mean(pos3(:,3)); %pos2(index2,3)
cpos4 = mean(pos4(:,3)); %pos2(index2,3)

tot_mean = (cpos1+cpos4)/2

cell = zeros(3);
cell(1,1:3) = str2num(data.textdata{8});
cell(2,1:3) = str2num(data.textdata{9});
cell(3,1:3) = str2num(data.textdata{10});
phi = acos(cell(1,1)/norm(cell(1,:)));
Rphi = [cos(-phi) -sin(-phi) 0;sin(-phi) cos(-phi) 0; 0 0 1];
pos1 = pos1*Rphi';
pos2 = pos2*Rphi';
pos3 = pos3*Rphi';
pos4 = pos4*Rphi';
cell = cell*Rphi';

pos1(:,3) = pos1(:,3) - tot_mean;
pos2(:,3) = pos2(:,3) - tot_mean;
pos3(:,3) = pos3(:,3) - tot_mean;
pos4(:,3) = pos4(:,3) - tot_mean;
cpos1 = mean(pos1(:,3));
cpos2 = mean(pos2(:,3));
cpos3 = mean(pos3(:,3));
cpos4 = mean(pos4(:,3));

for iat = 1 : nat1
    at_name1(iat) = string(data.textdata{10+iat});
end
for jat = 1 : nat2
    at_name2(jat) = string(data.textdata{10+nat1+jat});
end
for kat = 1 : nat3
    at_name3(kat) = string(data.textdata{10+nat1+nat2+kat});
end

for kat = 1 : nat4
    at_name4(kat) = string(data.textdata{10+nat1+nat2+nat3+kat});
end


ind1 = find(at_name1 == "Mo");
if(isempty(ind1))
    ind1 = find(at_name1 == "W");
end


ind2 = find(at_name2 == "Mo");
if(isempty(ind2))
    ind2 = find(at_name2 == "W");
end


ind3 = find(at_name3 == "Mo");
if(isempty(ind3))
    ind3 = find(at_name3 == "W");
end


ind4 = find(at_name4 == "Mo");
if(isempty(ind4))
    ind4 = find(at_name4 == "W");
end
% Create supercell for top layer
supercell = zeros(length(ind4)*36,3);
siat = 0;
for iiat = 1 : length(ind4)
    iat = ind4(iiat);
    for in = -3 : 2
        for jn = -3 : 2
            siat = siat + 1;
            supercell(siat,:) = pos4(iat,:) + in*cell(1,:) + jn*cell(2,:);
        end
    end
end
f4 = scatteredInterpolant(supercell(:,1),supercell(:,2),supercell(:,3),'natural');

supercell = zeros(length(ind3)*36,3);
siat = 0;
for iiat = 1 : length(ind3)
    iat = ind3(iiat);
    for in = -3 : 2
        for jn = -3 : 2
            siat = siat + 1;
            supercell(siat,:) = pos3(iat,:) + in*cell(1,:) + jn*cell(2,:);
        end
    end
end
f3 = scatteredInterpolant(supercell(:,1),supercell(:,2),supercell(:,3),'natural');

% Create supercell for middle layer
supercell = zeros(length(ind2)*36,3);
siat = 0;
for iiat = 1 : length(ind2)
    iat = ind2(iiat);
    for in = -3 : 2
        for jn = -3 : 2
            siat = siat + 1;
            supercell(siat,:) = pos2(iat,:) + in*cell(1,:) + jn*cell(2,:);
        end
    end
end
f2 = scatteredInterpolant(supercell(:,1),supercell(:,2),supercell(:,3),'natural');

supercell = zeros(length(ind1)*36,3);
siat = 0;
for iiat = 1 : length(ind1)
    iat = ind1(iiat);
    for in = -3 : 2
        for jn = -3 : 2
            siat = siat + 1;
            supercell(siat,:) = pos1(iat,:) + in*cell(1,:) + jn*cell(2,:);
        end
    end
end

f1 = scatteredInterpolant(supercell(:,1),supercell(:,2),supercell(:,3),'natural');

if(mode==1)
    min_z = zeros(length(ind1),1);
    for iiat = 1 : length(ind1)
        iat = ind1(iiat);
        dis = f2(pos1(iat,1),pos1(iat,2)) - f1(pos1(iat,1),pos1(iat,2));
        min_z(iiat) = dis;
    end
    
    out1 = min(min_z);
    out2 = max(min_z);

elseif(mode==2)
    if(class=="homos")
    %HOMO
    x2 = linspace(0.0,2*(cell(1,1)+cell(2,1)),500);
    y2 = linspace(0.0,2*(cell(1,2)+cell(2,2)),500);
    fullnorm = norm(cell(1,:)+cell(2,:));
    elseif(class=="heteros")
    % HETERO
    x2 = linspace(2*cell(2,1)-2*cell(1,1),0.0,200);
    y2 = linspace(2*cell(2,2),0.0,200);
    fullnorm = (norm(-cell(1,:)+cell(2,:)));
    end
    out1 = 0;
    out2 = 0;

[X2,Y2] = meshgrid(x2,y2);
Ztest1 = f1(x2,y2);
Ztest2 = f2(x2,y2);
Ztest3 = f3(x2,y2);
Ztest4 = f4(x2,y2);

l = sqrt(x2.^2+y2.^2);

plot(ax1,l/fullnorm,Ztest4,'Color',rgb,'Linewidth',2);
plot(ax2,l/fullnorm,Ztest3,'Color',rgb,'Linewidth',2);
plot(ax3,l/fullnorm,Ztest2,'Color',rgb,'Linewidth',2);
plot(ax4,l/fullnorm,Ztest1,'Color',rgb,'Linewidth',2);


elseif(mode==3)
supercell = zeros(length(ind1)*4,3);
siat = 0;
for iiat = 1 : length(ind1)
    iat = ind1(iiat);
    for in = 0 : 1
        for jn = 0 : 1
            siat = siat + 1;
            supercell(siat,:) = pos1(iat,:) + in*cell(1,:) + jn*cell(2,:);
        end
    end
end

f = figure;
set(f,'Position',[500.5 500.5 800 600]);
Zdiff = f3(supercell(:,1),supercell(:,2))-f2(supercell(:,1),supercell(:,2));
scatter(supercell(:,1),supercell(:,2),100,Zdiff,'filled');

hold on;
plot([0 cell(1,1)],[0 0],'k','LineWidth',2.0);
plot([0 cell(2,1)],[0 cell(2,2)],'k','LineWidth',2.0);
plot([cell(2,1) cell(2,1)+cell(1,1)],[cell(2,2) cell(2,2)],'k','LineWidth',2.0);
plot([cell(1,1) cell(1,1)+cell(2,1)],[0 cell(2,2)],'k','LineWidth',2.0);
plot([cell(1,1) 2*cell(1,1)],[0 0],'k','LineWidth',2.0);
plot([2*cell(1,1) 2*cell(1,1)+cell(2,1)],[0 cell(2,2)],'k','LineWidth',2.0);
plot([cell(1,1)+cell(2,1) 2*cell(1,1)+ cell(2,1) ],[cell(2,2) cell(2,2)],'k','LineWidth',2.0);
plot([cell(2,1) 2*cell(2,1)],[cell(2,2) 2*cell(2,2)],'k','LineWidth',2.0);
plot([2*cell(2,1) 2*cell(2,1)+cell(1,1)],[2*cell(2,2) 2*cell(2,2)],'k','LineWidth',2.0);
plot([2*cell(2,1)+cell(1,1) 2*cell(2,1)+2*cell(1,1)],[2*cell(2,2) 2*cell(2,2)],'k','LineWidth',2.0);
plot([cell(2,1)+cell(1,1) cell(2,1)+2*cell(1,1)],[cell(2,2) cell(2,2)],'k','LineWidth',2.0);
plot([cell(2,1)+cell(1,1) 2*cell(2,1)+cell(1,1)],[cell(2,2) 2*cell(2,2)],'k','LineWidth',2.0);
plot([cell(2,1)+2*cell(1,1) 2*cell(2,1)+2*cell(1,1)],[cell(2,2) 2*cell(2,2)],'k','LineWidth',2.0);
if(class=="homos")
% HOMO
 plot([0 2*cell(1,1)+2*cell(2,1)], [0 2*cell(2,2)],'g--','LineWidth',3.0);
 plot([55,105],[20,20],'w','LineWidth',5.0);
 text(70,25,'\textbf{50 \AA}','Interpreter','latex','Color','white','FontSize',28,'FontWeight','bold');
elseif(class=="heteros");
% HETERO
 plot([2*cell(1,1), 2*cell(2,1)], [0 2*cell(2,2)],'g--','LineWidth',3.0);
 plot([15,65],[10,10],'w','LineWidth',5.0);
 text(25,15,'\textbf{50 \AA}','Interpreter','latex','Color','white','FontSize',28);
end
axis([min(supercell(:,1))-25, max(supercell(:,1)) + 25, min(supercell(:,2))-50, max(supercell(:,2)) + 50]);
axis('off');
h = colorbar;ylabel(h,'ILS (\AA)','Interpreter','latex','FontSize',24);
colormap('hot');
x = get(gca,'position');
xh = get(h,'position');
xh(1)=0.15;
xh(2)=0.3;
xh(4)=0.5;
set(h,'position',xh);
set(gca,'position',x);
set(gca,'FontSize',28);

out1=0;
out2=0;
end
clear txt at_name1 at_name2 cpos1 cpos2 data natoms nat1 nat2 X Y tind1 ind1 ind2 tind1 ind1t tind2 ind3 data pos1 pos2 X2 Y2 x y x2 y2 f1 f2 Ztest1 Ztest2 supercell
