%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to set TB parameters
% Adapted from matlab Shiang Fang matlab code:
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017
% Reference citation: Phys. Rev. B 92, 205108 (2015).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onsite,type1,type2,type3,type4,type5,type6,lsocM,lsocX] = set_parameters(tmdc,gw_par)

% Read in parameters and classify them in terms of their type
TMDC_monolayer;
tbparms = mono_tbh_parm;
%vparms = mono_vint_parm;
%inter_par;

sq3=sqrt(3);
sq2=sqrt(2);

onsite = zeros(1,11);
type1 = zeros(11);
type2 = zeros(11);
type3 = zeros(11);
type4 = zeros(11);
type5 = zeros(11);
type6 = zeros(11);
% 40 Independent parameters

gw1 = 0.0;
gw2 = gw1;
if(gw_par)
   gw1 = 0.3624;
   gw2 = -0.2512;
end

% ON-SITE
ep1=gw1 + tbparms(1);
ep3=gw2 + tbparms(2);
ep4=gw2 + tbparms(3);
ep6=gw1 + tbparms(4);
ep7=gw1 + tbparms(5);
ep9=gw2 + tbparms(6);
ep10=gw2 + tbparms(7);
ep2=ep1;
ep5=ep4;
ep8=ep7;
ep11=ep10;
onsite = [ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8,ep9,ep10,ep11];

gw1 = 1.0;
gw2 = gw1;
if(gw_par)
   gw1 = 1.4209;
   gw2 = 1.1738;
end

% TYPE 1 OFF-DIAGONAL
type1(3,5)  = gw2*tbparms(19);
type1(6,8)  = gw1*tbparms(20);
type1(9,11) = gw2*tbparms(21);
type1(1,2)  = gw1*tbparms(22);
type1(3,4)  = gw2*tbparms(23);
type1(4,5)  = gw2*tbparms(24);
type1(6,7)  = gw1*tbparms(25);
type1(7,8)  = gw1*tbparms(26);
type1(9,10) = gw2*tbparms(27);
type1(10,11)= gw2*tbparms(28);

type1 = type1 + type1';

% TYPE 1 DIAGONAL
type1(1,1)  = gw1*tbparms(8);
type1(2,2)  = gw1*tbparms(9);
type1(3,3)  = gw2*tbparms(10);
type1(4,4)  = gw2*tbparms(11);
type1(5,5)  = gw2*tbparms(12);
type1(6,6)  = gw1*tbparms(13);
type1(7,7)  = gw1*tbparms(14);
type1(8,8)  = gw1*tbparms(15);
type1(9,9)  = gw2*tbparms(16);
type1(10,10)= gw2*tbparms(17);
type1(11,11)= gw2*tbparms(18);

gw1 = 1.0;
if(gw_par)
   gw1 = 1.0773;
end
% TYPE 5 OFF-DIAGONAL ONLY
type5(4,1) = gw1*tbparms(29);
type5(3,2) = gw1*tbparms(30);
type5(5,2) = gw1*tbparms(31);
type5(9,6) = gw1*tbparms(32);
type5(11,6)= gw1*tbparms(33);
type5(10,7)= gw1*tbparms(34);
type5(9,8)=  gw1*tbparms(35);
type5(11,8)= gw1*tbparms(36);

gw1 = 1.0;
if(gw_par)
   gw1 = 1.1871;
end
% TYPE 6 OFF-DIAGONAL ONLY
type6(9,6) = gw1*tbparms(37);
type6(11,6)= gw1*tbparms(38);
type6(9,8) = gw1*tbparms(39);
type6(11,8)= gw1*tbparms(40);


% generate all remaining parameters

% TYPE 2 OFF-DIAGONAL

type2(3,5)=(sq3/2)*type1(3,4)-0.5*type1(3,5);
type2(6,8)=(sq3/2)*type1(6,7)-0.5*type1(6,8);
type2(9,11)=(sq3/2)*type1(9,10)-0.5*type1(9,11);

type2(1,2)=(sq3/4)*(type1(1,1)-type1(2,2))-type1(1,2);
type2(4,5)=(sq3/4)*(type1(4,4)-type1(5,5))-type1(4,5);
type2(7,8)=(sq3/4)*(type1(7,7)-type1(8,8))-type1(7,8);
type2(10,11)=(sq3/4)*(type1(10,10)-type1(11,11))-type1(10,11);

type2(3,4)=0.5*type1(3,4)+(sq3/2)*type1(3,5);
type2(6,7)=0.5*type1(6,7)+(sq3/2)*type1(6,8);
type2(9,10)=0.5*type1(9,10)+(sq3/2)*type1(9,11);

type2 = type2+type2';

% TYPE 2 DIAGONAL
type2(1,1)=0.25*type1(1,1)+0.75*type1(2,2);
type2(4,4)=0.25*type1(4,4)+0.75*type1(5,5);
type2(7,7)=0.25*type1(7,7)+0.75*type1(8,8);
type2(10,10)=0.25*type1(10,10)+0.75*type1(11,11);

type2(2,2)=0.75*type1(1,1)+0.25*type1(2,2);
type2(5,5)=0.75*type1(4,4)+0.25*type1(5,5);
type2(8,8)=0.75*type1(7,7)+0.25*type1(8,8);
type2(11,11)=0.75*type1(10,10)+0.25*type1(11,11);

type2(3,3)=type1(3,3);
type2(6,6)=type1(6,6);
type2(9,9)=type1(9,9);

% TYPE 3 OFF-DIAGONAL ONLY

type3(3,5)=-(sq3/2)*type1(3,4)-0.5*type1(3,5);
type3(6,8)=-(sq3/2)*type1(6,7)-0.5*type1(6,8);
type3(9,11)=-(sq3/2)*type1(9,10)-0.5*type1(9,11);

type3(1,2)=-(sq3/4)*(type1(1,1)-type1(2,2))-type1(1,2);
type3(4,5)=-(sq3/4)*(type1(4,4)-type1(5,5))-type1(4,5);
type3(7,8)=-(sq3/4)*(type1(7,7)-type1(8,8))-type1(7,8);
type3(10,11)=-(sq3/4)*(type1(10,10)-type1(11,11))-type1(10,11);

type3(3,4)=0.5*type1(3,4)-(sq3/2)*type1(3,5);
type3(6,7)=0.5*type1(6,7)-(sq3/2)*type1(6,8);
type3(9,10)=0.5*type1(9,10)-(sq3/2)*type1(9,11);

% TYPE 4 OFF-DIAGONAL ONLY
type4(4,1)=0.25*type5(4,1)+0.75*type5(5,2);
type4(10,7)=0.25*type5(10,7)+0.75*type5(11,8);

type4(5,2)=0.75*type5(4,1)+0.25*type5(5,2);
type4(11,8)=0.75*type5(10,7)+0.25*type5(11,8);

type4(5,1)=-(sq3/4)*type5(4,1)+(sq3/4)*type5(5,2);
type4(4,2)=type4(5,1);
type4(11,7)=-(sq3/4)*type5(10,7)+(sq3/4)*type5(11,8);
type4(10,8)=type4(11,7);

type4(3,1)=-(sq3/2)*type5(3,2);
type4(9,7)=-(sq3/2)*type5(9,8);

type4(3,2)=-0.5*type5(3,2);
type4(9,8)=-0.5*type5(9,8);

type4(9,6)=type5(9,6);
type4(10,6)=-sq3*type5(11,6)/2;
type4(11,6)=-type5(11,6)/2;

type3 = type3+type3';
type4 = type4+type4';
type5 = type5+type5';
type6 = type6+type6';

%Vsigma=vparms(1);
%Vpi=vparms(2);
%R=vparms(3:4);
%eta=vparms(5:6);

lsocM = soc_m;
lsocX = soc_x;

end

