%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate k-point set for       %%
%% for bandstructure calculations %%
%% along the Gamma-M-K-Gamma path %%
%% and for DOS/LDOS/chern etc. on %%
%% a uniform M.P. grid            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                %% 
%% Written by Valerio Vitale      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_kpts,scale_axis,knum_tot,recL,frac_k] = generate_kpoints(task,multilayer,twisted,miniBZ,mb1,mb2,mb3,b1,b2,b3,knum,bse_shifted)

% If computing bandstructure set \Gamma-M-K-\Gamma path
if(task==1 || task==6)
    
    % If twisted bilayer, find high-symmetry k-points of monolayer BZ with respect to the vectors of the
    % mini-Brillouin zone
    if(multilayer && twisted && ~miniBZ)
       recL=[mb1;mb2;mb3];
       k1=[0,0,0];
       k2=b2/2;
       k3=(2*b1-b2)/3;
       %k1=(2*b1+b2)/3;
       %k2=(b1+2*b2)/6;
       %k3= [0 0 0];
       fM = k2*inv(recL);
       fK = k3*inv(recL);
       % Find closest M and K points
       C3 = [cos(pi/3) sin(pi/3) 0; -sin(pi/3) cos(pi/3) 0; 0 0 1];
       k3_all = [k3;k3*C3;k3*C3*C3;k3*C3*C3*C3;k3*C3*C3*C3*C3;k3*C3*C3*C3*C3*C3];
       [minv,iv] = min(sqrt(sum((k3_all - k2).^2,2)));
       k3 = k3_all(iv,:);
       fK = k3*inv(recL);
       if(abs(fM(1)) >= 1)
          if(fM(1)  > 0)
             fM(1) = fM(1) - floor(fM(1)); %
          else
             fM(1) = fM(1) - ceil(fM(1)); %
          end
       end
       if(abs(fM(2)) >= 1)
          if(fM(2)  > 0)
             fM(2) = fM(2) - floor(fM(2)); %
          else
             fM(2) = fM(2) - ceil(fM(2)); %
          end
       end
       if(abs(fK(1)) >= 1)
          if(fK(1)  > 0)
             fK(1) = fK(1) - floor(fK(1)); %
          else
             fK(1) = fK(1) - ceil(fK(1)); %
          end
       end
       if(abs(fK(2)) >= 1)
          if(fK(2)  > 0)
             fK(2) = fK(2) - floor(fK(2)); %
          else
             fK(2) = fK(2) - ceil(fK(2)); %
          end
       end
       fM = [0.5 0.0 0.0];
       fK = [2/3 -1/3 0];
       k2 = fM*recL;
       k3 = fK*recL;
    else
       recL = [mb1;mb2;mb3];%
       k1 = [0,0,0];
       k2 = [0.5,0,0]*recL;
       k3 = [2/3,-1/3,0]*recL;
    end
    % Generate path in k-space
    
    scan_klist=[k1;k2;k3;k1];
    
    [ all_kpts, scale_axis] = generate_k_line( knum, scan_klist );
    dk = norm(all_kpts(2,:) - all_kpts(1,:));
    dk = dk*0.529177210903;
    knum_tot=size(all_kpts);
    knum_tot=knum_tot(1);
    frac_k = all_kpts*(inv(recL));
    %[all_kpts, ones(knum_tot,1)/knum_tot]
elseif(task==2 || task==3 || task==4 || task==5 || task==7)
    
    recL = [mb1;mb2;mb3];%
    K = zeros(knum*knum,3);
    all_kpts = zeros(knum*knum,3);
    kxi=0;
    kyi=0;
    kxf=(knum-1)/knum;
    kyf=kxf;
    kx=linspace(kxi,kxf,knum);
    ky=linspace(kyi,kyf,knum);
    ikk = 0;
    for ikx = 1 : knum
       for iky = 1 : knum
           ikk = ikk + 1;
           K(ikk,:) = [kx(ikx),ky(iky),0.0];
       end
    end
    scale_axis = sqrt(sum(K.^2,2));
    all_kpts = K*recL;
    knum_tot = knum*knum;
    % Shift the uniform grid by a small dk with a positive random direction
    % for the calculation of BSE eigenvalues/eigenfunctions
    if(task == 5 && bse_shifted)
       % Fix the seed for the random number generator
       rng(42);
       cost = rand(1);
       sint = sqrt(1-cost);
       shift = norm(all_kpts(2,:)-all_kpts(1,:))/10 * [cost,sint,0];
       all_kpts = all_kpts + shift;
    end
    clear K;
end
end
