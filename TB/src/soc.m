%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to compute SOC matrix in 11-orbital 
% basis:
% <\phi_{i,sigma} | \lambda L.S | \phi_{j,sigma'}>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% lsocM - SO coupling strength for metal atoms
% lsocX - SO coupling strength for chalcogen atoms
%
% Outputs:
% Hsoc - 22x22 SOC matrix in 22-orbital basis
% cjmd - Clebsh-Gordan unitary matrix for d-orbitals
% cjmp - Clebsh-Gordan unitary matrix for p-orbitals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Valerio Vitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hsoc,Up,Ud]= soc(lsocM,lsocX)
% d-orbitals first
l = 2;
s = 1/2;
ml = [-l, -(l-1), 0, (l-1), l];
ms = [-s, s];
jmax = l+s;
jmin = l-s;
Mjmax = [jmax, jmax-1, jmax-2, -(jmax-2), -(jmax-1), -jmax]; % -J, ..., J
Mjmin = [jmin, jmin-1, -(jmin-1), -jmin]; % -J, ..., J
M = [size(Mjmax,2),size(Mjmin,2)];
MM = [Mjmax, Mjmin];
J = [jmax,jmin]; % l+s,..,|l-s|
SOC_M = zeros(10);

% Generates SOC matrix in |JMls> basis
i = 0;
for j = 1 : 2
    lj = J(j);
    for m = 1 : M(j)
        i = i + 1;
        SOC_M(i,i) = 0.5*lsocM*(lj*(lj+1) - l*(l+1) -s*(s+1));
    end
end
% Clebsh-Gordan
p = 0;
cjmd = zeros(10);
for n = 1 : size(ms,2)
    lms = ms(n);
    for m = 1 : size(ml,2)
        lml = ml(m);
        p = p + 1;
        q = 0;
        for j = 1 : 2
            lJ = J(j);
            for mm = 1 : M(j)
                q = q + 1;
                lM = MM(mm + (j-1)*M(1));
                cjmd(p,q) = clebsh_gordan_d(lJ,lM,lml,lms);
            end
        end
    end
end

% then p-orbitals
l = 1;
s = 1/2;
ml = [-l, 0, l];
ms = [-s, s];
jmax = l+s;
jmin = l-s;
Mjmax = [jmax, jmax-1, -(jmax-1), -jmax]; % -J, ..., J
Mjmin = [jmin, -jmin]; % -J, ..., J
M = [size(Mjmax,2),size(Mjmin,2)];
MM = [Mjmax, Mjmin];
J = [jmax,jmin]; % l+s,..,|l-s|
SOC_X = zeros(6);

% Generates SOC matrix in |JMls> basis
i = 0;
for j = 1 : 2
    lj = J(j);
    for m = 1 : M(j)
        i = i + 1;
        SOC_X(i,i) = 0.5*lsocX*(lj*(lj+1) - l*(l+1) -s*(s+1));
    end
end
% Clebsh-Gordan
p = 0;
cjmp = zeros(6);
for n = 1 : size(ms,2)
    lms = ms(n);
    for m = 1 : size(ml,2)
        lml = ml(m);
        p = p + 1;
        q = 0;
        for j = 1 : 2
            lJ = J(j);
            for mm = 1 : M(j)
                q = q + 1;
                lM = MM(mm + (j-1)*M(1));
                cjmp(p,q) = clebsh_gordan_p(lJ,lM,lml,lms);
            end
        end
    end
end

% Convert from spherical harmonics to cubic harmonics
sq2 = sqrt(2);
ang_1 = [0,1,0,0,0,0; %pz_down
    -1/sq2,0,1/sq2,0,0,0; %px_down
    1i/sq2,0,1i/sq2,0,0,0;
    0,0,0,0,1,0; %pz_up
    0,0,0,-1/sq2,0,1/sq2; %px_up
    0,0,0,1i/sq2,0,1i/sq2]; %py_up


ang_2 = [0,-1/sq2,0,1/sq2,0,0,0,0,0,0;
    0,1i/sq2,0,1i/sq2,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    -1i/sq2,0,0,0,1i/sq2,0,0,0,0,0;
    1/sq2,0,0,0,1/sq2,0,0,0,0,0;
    0,0,0,0,0,0,-1/sq2,0,1/sq2,0; %dxz up
    0,0,0,0,0,0,1i/sq2,0,1i/sq2,0; %dyz up
    0,0,0,0,0,0,0,1,0,0; %%dz2 up
    0,0,0,0,0,-1i/sq2,0,0,0,1i/sq2; %dxy up
    0,0,0,0,0,1/sq2,0,0,0,1/sq2]; %dx2-y2 up
% Evaluate SOC matrices from Clebsh-Gordan and spherical to cubic matrices
Up = cjmp'*(conj(ang_1'));
Ud = cjmd'*(conj(ang_2'));
%Up = inv(cjmp)*inv(ang_1);
%Ud = inv(cjmd)*inv(ang_2);
Hsocp = Up'*SOC_X*Up;
Hsocd = Ud'*SOC_M*Ud;

Hsoc = zeros(22);

% Deposit SOC matrices into Hsoc matrix
Hsoc(1:2,1:2)     = Hsocd(1:2,1:2);
Hsoc(1:2,6:8)     = Hsocd(1:2,3:5);
Hsoc(1:2,12:13)   = Hsocd(1:2,6:7);
Hsoc(1:2,17:19)   = Hsocd(1:2,8:10);
Hsoc(6:8,1:2)     = Hsocd(3:5,1:2);
Hsoc(6:8,6:8)     = Hsocd(3:5,3:5);
Hsoc(6:8,12:13)   = Hsocd(3:5,6:7);
Hsoc(6:8,17:19)   = Hsocd(3:5,8:10);
Hsoc(12:13,1:2)   = Hsocd(6:7,1:2);
Hsoc(12:13,6:8)   = Hsocd(6:7,3:5);
Hsoc(12:13,12:13) = Hsocd(6:7,6:7);
Hsoc(12:13,17:19) = Hsocd(6:7,8:10);
Hsoc(17:19,1:2)   = Hsocd(8:10,1:2);
Hsoc(17:19,6:8)   = Hsocd(8:10,3:5);
Hsoc(17:19,12:13) = Hsocd(8:10,6:7);
Hsoc(17:19,17:19) = Hsocd(8:10,8:10);

Hsoc(3:5,3:5)    = Hsocp(1:3,1:3);
Hsoc(3:5,14:16)  = Hsocp(1:3,4:6);
Hsoc(14:16,3:5)  = Hsocp(4:6,1:3);
Hsoc(14:16,14:16) = Hsocp(4:6,4:6);

Hsoc(9:11,9:11)  = Hsocp(1:3,1:3);
Hsoc(9:11,20:22) = Hsocp(1:3,4:6);
Hsoc(20:22,9:11) = Hsocp(4:6,1:3);
Hsoc(20:22,20:22) = Hsocp(4:6,4:6);

%new basis:

% noddxu=4;
% noddxd=15;
% noddyu=5;
% noddyd=16;
% noddzu=3;
% noddzd=14;
% noddxzu=1;
% noddxzd=12;
% noddyzu=2;
% noddyzd=13;
% nevenxu=10;
% nevenxd=21;
% nevenyu=11;
% nevenyd=22;
% nevenzu=9;
% nevenzd=20;
% nevenz2u=6;
% nevenz2d=17;
% nevenxyu=7;
% nevenxyd=18;
% nevenx2y2u=8;
% nevenx2y2d=19;
% 
% %old basis:
% odz2u=1;
% odz2d=12;
% odxyu=2;
% odxyd=13;
% odx2y2u=3;
% odx2y2d=14;
% odxzu=4;
% odxzd=15;
% odyzu=5;
% odyzd=16;
% os1xu=6;
% os1xd=17;
% os1yu=7;
% os1yd=18;
% os1zu=8;
% os1zd=19;
% os2xu=9;
% os2xd=20;
% os2yu=10;
% os2yd=21;
% os2zu=11;
% os2zd=22;
% 
% mat_wan_oe=zeros(22,22);
% 
% mat_wan_oe(odxzu,noddxzu)=1.0;
% mat_wan_oe(odxzd,noddxzd)=1.0;
% mat_wan_oe(odyzu,noddyzu)=1.0;
% mat_wan_oe(odyzd,noddyzd)=1.0;
% mat_wan_oe(odz2u,nevenz2u)=1.0;
% mat_wan_oe(odz2d,nevenz2d)=1.0;
% mat_wan_oe(odxyu,nevenxyu)=1.0;
% mat_wan_oe(odxyd,nevenxyd)=1.0;
% mat_wan_oe(odx2y2u,nevenx2y2u)=1.0;
% mat_wan_oe(odx2y2d,nevenx2y2d)=1.0;
% 
% mat_wan_oe(os1xu,noddxu)=1/sq2;
% mat_wan_oe(os1xd,noddxd)=1/sq2;
% mat_wan_oe(os2xu,noddxu)=-1/sq2;
% mat_wan_oe(os2xd,noddxd)=-1/sq2;
% mat_wan_oe(os1yu,noddyu)=1/sq2;
% mat_wan_oe(os1yd,noddyd)=1/sq2;
% mat_wan_oe(os2yu,noddyu)=-1/sq2;
% mat_wan_oe(os2yd,noddyd)=-1/sq2;
% mat_wan_oe(os1zu,noddzu)=1/sq2;
% mat_wan_oe(os1zd,noddzd)=1/sq2;
% mat_wan_oe(os2zu,noddzu)=1/sq2;
% mat_wan_oe(os2zd,noddzd)=1/sq2;
% 
% mat_wan_oe(os1xu,nevenxu)=1/sq2;
% mat_wan_oe(os1xd,nevenxd)=1/sq2;
% mat_wan_oe(os2xu,nevenxu)=1/sq2;
% mat_wan_oe(os2xd,nevenxd)=1/sq2;
% mat_wan_oe(os1yu,nevenyu)=1/sq2;
% mat_wan_oe(os1yd,nevenyd)=1/sq2;
% mat_wan_oe(os2yu,nevenyu)=1/sq2;
% mat_wan_oe(os2yd,nevenyd)=1/sq2;
% mat_wan_oe(os1zu,nevenzu)=1/sq2;
% mat_wan_oe(os1zd,nevenzd)=1/sq2;
% mat_wan_oe(os2zu,nevenzu)=-1/sq2;
% mat_wan_oe(os2zd,nevenzd)=-1/sq2;

%Hsoc = mat_wan_oe*Hsoc*mat_wan_oe';
end

function CJM = clebsh_gordan_d(J,M,ml,ms)
  CJM = 0.0;
  if (J==5/2 && M == 5/2)
     if (ml == 2 && ms == 1/2)
        CJM = 1;
     end
  elseif (J==5/2 && M == 3/2)
     if (ml == 2 && ms == -1/2)
        CJM = sqrt(1/5);
     elseif( ml == 1 && ms == 1/2)
        CJM = sqrt(4/5);
     end
  elseif( J==5/2 && M == 1/2)
     if (ml == 1 && ms == -1/2)
        CJM = sqrt(2/5);
     elseif (ml == 0 && ms == 1/2)
        CJM = sqrt(3/5);
     end
  elseif( J==5/2 && M == -1/2)
     if (ml == 0 && ms == -1/2)
        CJM = sqrt(3/5);
     elseif (ml ==-1 && ms == 1/2)
        CJM = sqrt(2/5);
     end
  elseif( J==5/2 && M == -3/2)
     if (ml ==-1 && ms == -1/2)
        CJM = sqrt(4/5);
     elseif (ml ==-2 && ms == 1/2)
        CJM = sqrt(1/5);
     end
  elseif (J==5/2 && M == -5/2)
     if (ml ==-2 && ms ==-1/2)
        CJM = 1;
     end
  elseif (J==3/2 && M == 3/2)
     if (ml == 2 && ms == -1/2)
        CJM = sqrt(4/5);
     elseif (ml == 1 && ms == 1/2)
        CJM = -sqrt(1/5);
     end
  elseif (J==3/2 && M == 1/2)
     if (ml == 1 && ms == -1/2)
        CJM = sqrt(3/5);
     elseif (ml == 0 && ms == 1/2)
        CJM = -sqrt(2/5);
     end
  elseif (J==3/2 && M ==-1/2)
     if (ml == 0 && ms == -1/2)
        CJM = sqrt(2/5);
     elseif (ml ==-1 && ms == 1/2)
        CJM = -sqrt(3/5);
     end
  elseif (J==3/2 && M ==-3/2)
     if (ml ==-1 && ms == -1/2)
        CJM = sqrt(1/5);
     elseif (ml ==-2 && ms == 1/2)
        CJM = -sqrt(4/5);
     end
  else
     CJM = 0.0;
  end
end


function CJM = clebsh_gordan_p(J,M,ml,ms)
  CJM = 0.0;
  sq2=sqrt(2);
  sq3=sqrt(3);
  if (J==3/2 && M == 3/2)
     if (ml == 1 && ms == 1/2)
        CJM = 1;
     end
  elseif (J==3/2 && M == 1/2)
     if (ml == 1 && ms == -1/2)
        CJM = 1/sq3;
     elseif (ml == 0 && ms == 1/2)
        CJM = sq2/sq3;
     end
  elseif (J==1/2 && M ==1/2)
     if (ml == 1 && ms == -1/2)
        CJM = sq2/sq3;
     elseif (ml ==0 && ms == 1/2)
        CJM = -1/sq3;
     end
  elseif (J==3/2 && M ==-1/2)
     if (ml ==0 && ms == -1/2)
        CJM = sq2/sq3;
     elseif (ml ==-1 && ms == 1/2)
        CJM = 1/sq3;
     end
  elseif (J==1/2 && M==-1/2)
     if (ml==0 && ms == -1/2)
         CJM = 1/sq3;
     elseif(ml == -1 && ms == 1/2)
         CJM = -sq2/sq3;
     end
  elseif (J==3/2 && M == -3/2)
      if (ml == -1 && ms == -1/2)
          CJM = 1;
      end
  else
     CJM = 0.0;
  end
end
