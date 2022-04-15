function write_kpts(npts,ralat,recL,th_string,str_m,str_n,path,bilayer)

        if(bilayer==1)
           Gamma = 0*ralat(:,1) + 0*ralat(:,2);       
           K1 = 2/3*ralat(:,1) + 1/3*ralat(:,2);
           K2 = 1/3*ralat(:,1) + 2/3*ralat(:,2);
           M = 1/2*ralat(:,1) + 1/2*ralat(:,2);
        elseif(bilayer==2)
           Gamma = 0*ralat(:,1) + 0*ralat(:,2);       
           K1 = 2/3*ralat(:,1) - 1/3*ralat(:,2);
           K2 = -2/3*ralat(:,1) + 1/3*ralat(:,2);
           M = 0.5*ralat(:,1);
        end
        
        % Generate kpoints.in along Gamma-K-M-Gamma path or K-Gamma-M-K'
        kpt_name = "kpoints";
        kpt_name = join([kpt_name,th_string],"_theta=");
        kpt_name = join([kpt_name,str_m],"_m=");
        kpt_name = join([kpt_name,str_n],"_n=");
        kpath_name = "kpaths";
        kpath_name = join([kpath_name,th_string],"_theta=");
        kpath_name = join([kpath_name,str_m],"_m=");
        kpath_name = join([kpath_name,str_n],"_n=");
        fileID = fopen(kpt_name,'w');
        fileID2= fopen(kpath_name,'w');
        
        disp('    Selected path in k-space:')
        tpath = join(['    ',path]);
        disp(tpath)
        switch path
            case 'KGMKp'
                
                disp("    K-Gamma-M-K'")
                printkpts_KMGK(npts,K1,K2,Gamma,M,recL,fileID,fileID2);
            case 'MGKM'
                
                disp('    M-Gamma-K-M')
                printkpts_MGKM(npts,K1,Gamma,M,recL,fileID,fileID2);
            case 'GMKG'
                
                disp('    Gamma-M-K-Gamma')
                printkpts_GMKG(npts,K1,Gamma,M,recL,fileID,fileID2);
            case 'GKM'
                
                disp('    Gamma-K-M')
                printkpts_GKM(npts,K1,Gamma,M,recL,fileID,fileID2);
            case 'KM'
                
                disp('    K-M')
                printkpts_KM(npts,K1,M,recL,fileID,fileID2);
            case 'uniform'
                
                disp('    uniform MP grid')
            otherwise
                error('Path in k-space not supported. Available paths are: K-Gamma-M-Kprime (KGMKp), M-Gamma-K-M (MGKM), Gamma-M-K-Gamma (GMKG), Gamma-K-M (GKM), K-M (KM) and Monkhorst-Pack (uniform)')
        end
        
        fclose(fileID);
        fclose(fileID2);
end

function printkpts_GKM(npts,K,Gamma,M,recL,fileID,fileID2)

        
        fGamma = Gamma'*inv(recL);
        fK = K'*inv(recL);
        fM = M'*inv(recL);
        fGamma = fGamma - floor(fGamma);
        fK = frac3(fK);
        fM = frac3(fM);
        distGK = sqrt((K(1)-Gamma(1))^2 + (K(2)-Gamma(2))^2);
        distMK = sqrt((M(1)-K(1))^2 + (M(2)-K(2))^2);
        
        % Find total number of K points
        npts1 = floor(distGK*npts/distMK);
        nkpts = npts + npts1;
        
        fprintf(fileID,"%i\n",nkpts);
        fprintf(fileID2,"%2.6f %2.6f\n", M(1), M(2));
        fprintf(fileID2,"%2.6f %2.6f\n", K(1), K(2));
        dk = abs(distGK/npts);                                              
        ki = [fGamma(1),fGamma(2)];
        kf = [fK(1),fK(2)];
        dist_tot = 0;
        kx = linspace(ki(1),kf(1),npts);
        ky = linspace(ki(2),kf(2),npts);
       
        for ik = 1 : npts
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik), ky(ik), 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n", dist_tot + ik*dk);
        end
        % K -> K'
        dk = abs(distMK/npts1);
        ki = [fK(1),fK(2)];
        kf = [fM(1),fM(2)];
        kx = linspace(ki(1),kf(1),npts1+1);
        ky = linspace(ki(2),kf(2),npts1+1);
        dist_tot = dist_tot + distGK;
        for ik = 2 : npts1+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end

end

function printkpts_KM(npts,K,M,recL,fileID,fileID2)
        fK = K'*inv(recL);
        fM = M'*inv(recL);
        fK = frac3(fK);
        fM = frac3(fM);
        distMK = sqrt((M(1)-K(1))^2 + (M(2)-K(2))^2);
        
        % Find total number of K points
        nkpts = npts;
        fprintf(fileID,"%i\n",nkpts);
        fprintf(fileID2,"%2.6f %2.6f\n", M(1), M(2));
        fprintf(fileID2,"%2.6f %2.6f\n", K(1), K(2));
        dk = abs(distMK/npts);
        ki = [fK(1),fK(2)];
        kf = [fM(1),fM(2)];
        kx = linspace(ki(1),kf(1),npts);
        ky = linspace(ki(2),kf(2),npts);
        dist_tot = 0.0;
        for ik = 1 : npts
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + ik*dk);
        end
end

function printkpts_KMGK(npts,K1,K2,Gamma,M,recL,fileID,fileID2)      
        fGamma = Gamma'*inv(recL);
        fK1 = K1'*inv(recL);
        fK2 = K2'*inv(recL);
        fM = M'*inv(recL);
        fGamma = fGamma - floor(fGamma);
        fK1 = frac3(fK1);
        fK2 = frac3(fK2);
        fM = frac3(fM);
        distK2G = sqrt((K2(1)-Gamma(1))^2 + (K2(2)-Gamma(2))^2);
        distGM = sqrt((M(1)-Gamma(1))^2 + (M(2)-Gamma(2))^2);
        distMK = sqrt((M(1)-K1(1))^2 + (M(2)-K1(2))^2);
        
        % Find total number of K points
        npts1 = floor(distGM*npts/distK2G);
        npts2 = floor(distMK*npts/distK2G);
        nkpts = npts + npts1 + npts2;
        
        fprintf(fileID,"%i\n",nkpts);
        fprintf(fileID2,"%2.6f %2.6f\n", K2(1), K2(2));
        fprintf(fileID2,"%2.6f %2.6f\n", M(1), M(2));
        fprintf(fileID2,"%2.6f %2.6f\n", K1(1), K1(2));
        dk = abs(distK2G/npts);                                              
        ki = [fK2(1),fK2(2)];
        kf = [fGamma(1),fGamma(2)];
        dist_tot = 0;
        kx = linspace(ki(1),kf(1),npts);
        ky = linspace(ki(2),kf(2),npts);
       
        for ik = 1 : npts
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik), ky(ik), 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n", dist_tot + ik*dk);
        end
        % K -> K'
        dk = abs(distGM/npts1);
        ki = [fGamma(1),fGamma(2)];
        kf = [fM(1),fM(2)];
        kx = linspace(ki(1),kf(1),npts1+1);
        ky = linspace(ki(2),kf(2),npts1+1);
        dist_tot = dist_tot + distK2G;
        for ik = 2 : npts1+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
        % K' -> Gamma
        dk = abs(distMK/npts2);
        ki = [fM(1),fM(2)];
        kf = [fK1(1),fK1(2)];
        kx = linspace(ki(1),kf(1),npts2+1);
        ky = linspace(ki(2),kf(2),npts2+1);
        dist_tot = dist_tot + distGM;
        for ik = 2 : npts2+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
end

function printkpts_MGKM(npts,K1,Gamma,M,recL,fileID,fileID2)      
        fGamma = Gamma'*inv(recL);
        fK1 = K1'*inv(recL);
        fM = M'*inv(recL);
        fGamma = fGamma - floor(fGamma);
        fK1 = frac3(fK1);
        fM = frac3(fM);
        distGM = sqrt((M(1)-Gamma(1))^2 + (M(2)-Gamma(2))^2);
        distKG = sqrt((K1(1)-Gamma(1))^2 + (K1(2)-Gamma(2))^2);
        distMK = sqrt((M(1)-K1(1))^2 + (M(2)-K1(2))^2);
        
        % Find total number of K points
        npts1 = floor(distKG*npts/distGM);
        npts2 = floor(distMK*npts/distGM);
        nkpts = npts + npts1 + npts2;
        
        fprintf(fileID,"%i\n",nkpts);
        fprintf(fileID2,"%2.6f %2.6f\n", M(1), M(2));
        fprintf(fileID2,"%2.6f %2.6f\n", K1(1), K1(2));
        dk = abs(distGM/npts);                                              
        ki = [fM(1),fM(2)];
        kf = [fGamma(1),fGamma(2)];
        dist_tot = 0;
        kx = linspace(ki(1),kf(1),npts);
        ky = linspace(ki(2),kf(2),npts);
       
        for ik = 1 : npts
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik), ky(ik), 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n", dist_tot + ik*dk);
        end
        % K -> K'
        dk = abs(distKG/npts1);
        ki = [fGamma(1),fGamma(2)];
        kf = [fK1(1),fK1(2)];
        kx = linspace(ki(1),kf(1),npts1+1);
        ky = linspace(ki(2),kf(2),npts1+1);
        dist_tot = dist_tot + distGM;
        for ik = 2 : npts1+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
        % K' -> Gamma
        dk = abs(distMK/npts2);
        ki = [fK1(1),fK1(2)];
        kf = [fM(1),fM(2)];
        kx = linspace(ki(1),kf(1),npts2+1);
        ky = linspace(ki(2),kf(2),npts2+1);
        dist_tot = dist_tot + distKG;
        for ik = 2 : npts2+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
end

function printkpts_GMKG(npts,K1,Gamma,M,recL,fileID,fileID2)      
        fGamma = Gamma'*inv(recL);
        fK1 = K1'*inv(recL);
        fM = M'*inv(recL);
        fGamma = fGamma - floor(fGamma);
        fK1 = frac3(fK1)
        fM = frac3(fM)
        distGM = sqrt((M(1)-Gamma(1))^2 + (M(2)-Gamma(2))^2);
        distMK = sqrt((M(1)-K1(1))^2 + (M(2)-K1(2))^2);
        distKG = sqrt((K1(1)-Gamma(1))^2 + (K1(2)-Gamma(2))^2);
        
        % Find total number of K points
        npts1 = floor(distMK*npts/distGM);
        npts2 = floor(distKG*npts/distGM);
        nkpts = npts + npts1 + npts2;
        
        fprintf(fileID,"%i\n",nkpts);
        fprintf(fileID2,"%2.6f %2.6f\n", M(1), M(2));
        fprintf(fileID2,"%2.6f %2.6f\n", K1(1), K1(2));
        dk = abs(distGM/npts)                                              
        ki = [fGamma(1),fGamma(2)];
        kf = [fM(1),fM(2)];
        dist_tot = 0;
        kx = linspace(ki(1),kf(1),npts);
        ky = linspace(ki(2),kf(2),npts);
       
        for ik = 1 : npts
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik), ky(ik), 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n", dist_tot + ik*dk);
        end
        % K -> K'
        dk = abs(distMK/npts1);
        ki = [fM(1),fM(2)];
        kf = [fK1(1),fK1(2)];
        kx = linspace(ki(1),kf(1),npts1+1);
        ky = linspace(ki(2),kf(2),npts1+1);
        dist_tot = dist_tot + distGM;
        for ik = 2 : npts1+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
        % K' -> Gamma
        dk = abs(distKG/npts2);
        ki = [fK1(1),fK1(2)];
        kf = [fGamma(1),fGamma(2)];
        kx = linspace(ki(1),kf(1),npts2+1);
        ky = linspace(ki(2),kf(2),npts2+1);
        dist_tot = dist_tot + distMK;
        for ik = 2 : npts2+1
           fprintf(fileID,"%2.6f   %2.6f   %2.6f\n", kx(ik) , ky(ik) , 1.0/nkpts);
           fprintf(fileID2,"%2.6f\n",dist_tot + (ik-1)*dk);
        end
end

function y = frac(fK)
        if (fK(1) > 0)
            y(1) = fK(1) - round(fK(1)-0.5);
        elseif(fK(1) < 0)
            y(1) = fK(1) - round(fK(1)+0.5);
        elseif(fK(1) == 0)
            y(1) = 0;
        end
        if (fK(2) > 0)
            y(2) = fK(2) - round(fK(2)-0.5);
        elseif (fK(2) < 0)
            y(2) = fK(2) - round(fK(2)+0.5);
        elseif(fK(2) == 0)
            y(2) = 0;
        end
end

function y = frac3(fK)
        a = round(frac(fK),3);
        if (abs(a(1)) > 1/2 && abs(a(2)) >= 1/2)
            if (a(1) > 0 && a(2) > 0)
                y(1) = a(1);
                y(2) = round(a(2) - 1.0,3);
            elseif(a(1) < 0 && a(2) < 0)
                y(1) = round(a(1) + 1.0,3);
                y(2) = a(2);
            else
                y=a;
            end
        else
            y=a;
        end
end
