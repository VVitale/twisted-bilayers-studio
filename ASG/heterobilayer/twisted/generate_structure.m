function [nat,nat1,nat2,pos,pos2,id1,id2,at_num1,at_num2] = generate_structure(a,u,c1,c2,z,I,Rtheta,ncell,alat1,ialatL,tol,pos,pos2,id1,id2,id11,id12,id13,id21,id22,id23,at_num1,at_num2,at_num11,at_num12,at_num13,at_num21,at_num22,at_num23,layer1,layer2,phase)
        k1=0;
        if(layer1 == "TMD" )
           for nj = -ncell:ncell
               for ni = -ncell:ncell
                   r1 = ni*alat1(:,1) + nj*alat1(:,2);% + [0.2;0.2]; 
                   r2 = r1 + 1/3*(2*alat1(:,1)+alat1(:,2));%[ r1(1); r1(2)+a/sqrt(3)];
                   r3 = r2;
                   l1 = ialatL*r1;
                   l2 = ialatL*r2;
                   % Check if atom coordinates is in moire cell
                   if (l1(1) < (1-tol) && l1(1) >= -tol && l1(2) < (1-tol) && l1(2) >= -tol) 
                       k1 = k1 + 1;
                       pos(k1,1)=r1(1);
                       pos(k1,2)=r1(2);
                       pos(k1,3)=0.25*c1;
                       id1(k1) = id11;
                       at_num1(k1) = at_num11;
                       lprint1(k1) = true;
                   end
                   if (l2(1) < (1-tol) && l2(1) >= -tol && l2(2) < (1-tol) && l2(2) >= -tol)
                       k1 = k1 + 1;
                       pos(k1,1)=r2(1);
                       pos(k1,2)=r2(2);
                       pos(k1,3)= -(0.5-u)*c1;
                       id1(k1) = id12;
                       at_num1(k1) = at_num12;
                       lprint1(k1) = true;
                       k1 = k1 + 1;
                       pos(k1,1)=r3(1);
                       pos(k1,2)=r3(2);
                       pos(k1,3)= (1-u)*c1;
                       id1(k1) = id13;
                       at_num1(k1) = at_num13;
                       lprint1(k1) = true;
                   end
               end
           end
        elseif(layer1 == "graphene")
           for nj = -ncell:ncell
               for ni = -ncell:ncell
                   r1 = ni*alat1(:,1) + nj*alat1(:,2); 
                   r2 = [ r1(1); r1(2)+a/sqrt(3)];
                   l1 = ialatL*r1;
                   l2 = ialatL*r2;
                   % Check if atom coordinates is in moire cell
                   if (l1(1) < (1-tol) && l1(1) >= -tol && l1(2) < (1-tol) && l1(2) >= -tol) 
                       k1 = k1 + 1;
                       pos(k1,1)=r1(1);
                       pos(k1,2)=r1(2);
                       pos(k1,3)=0.0;
                       id1(k1) = id11;
                       at_num1(k1) = at_num11;
                       lprint1(k1) = true;
                   end
                   if (l2(1) < (1-tol) && l2(1) >= -tol && l2(2) < (1-tol) && l2(2) >= -tol)
                       k1 = k1 + 1;
                       pos(k1,1)=r2(1);
                       pos(k1,2)=r2(2);
                       pos(k1,3)=0.0;
                       id1(k1) = id12;
                       at_num1(k1) = at_num12;
                       lprint1(k1) = true;
                   end
               end
           end
        end
        nat1 = k1;

        % Top layer
        k2 = 0;
        if(layer2 == "TMD")
           for ni=-ncell:ncell
               for nj = -ncell:ncell
                   r1 = ni*alat1(:,1) + nj*alat1(:,2);% + [0.2; 0.2];
                   r2 = r1 + 1/3*(2*alat1(:,1)+alat1(:,2));% [r1(1); r1(2)+a/sqrt(3)];
                   if(phase=='AA')
                      rr1 = I*Rtheta*r1;
                      rr2 = I*Rtheta*r2;
                      rr3 = rr2;
                   elseif(phase=='AB')
                      rr1 = I*Rtheta*(-r1 + [0.0; a]);%/sqrt(3)]);
                      rr2 = I*Rtheta*(-r2 + [0.0; a]);%/sqrt(3)]);
                      rr3 = rr2;
                   end

                   l1 = ialatL*rr1;
                   l2 = ialatL*rr2;
                   % Check atom belongs to moirè cell within small
                   % tollerance
                   if (l1(1) < (1-tol) && l1(1) >= -tol && l1(2) < (1-tol) && l1(2) >= -tol) 
                       k2 = k2 + 1;
                       pos2(k2,1) = rr1(1);
                       pos2(k2,2) = rr1(2);
                       pos2(k2,3) = 0.75*c2;
                       id2(k2) = id21;
                       at_num2(k2) = at_num21;
                       lprint2(k2) = true;
                       
                   end
                   if (l2(1) < (1-tol) && l2(1) >= -tol && l2(2) < (1-tol) && l2(2) >= -tol)
                       k2 = k2 + 1;
                       pos2(k2,1) = rr2(1);
                       pos2(k2,2) = rr2(2);
                       pos2(k2,3) = u*c2;
                       id2(k2) = id22;
                       at_num2(k2) = at_num22;
                       lprint2(k2) = true;
                       k2 = k2 + 1;
                       pos2(k2,1) = rr3(1);
                       pos2(k2,2) = rr3(2);
                       pos2(k2,3) = (1.5-u)*c2;
                       id2(k2) = id23;
                       at_num2(k2) = at_num23;
                       lprint2(k2) = true;
                   end
               end
           end
        elseif(layer2 == "graphene")
           for ni=-ncell:ncell
               for nj = -ncell:ncell
                   r1 = ni*alat1(:,1) + nj*alat1(:,2);
                   r2 = [r1(1); r1(2)+a/sqrt(3)];
                   if(phase=='AA')
                      rr1 = I*Rtheta*r1;
                      rr2 = I*Rtheta*r2;
                   elseif(phase=='AB')
                      rr1 = I*Rtheta*r2;
                      rr2 = I*Rtheta*[r2(1);r(2)+a/sqrt(3)];
                   end
                   l1 = ialatL*rr1;
                   l2 = ialatL*rr2;
                   % Check atom belongs to moirè cell within small
                   % tollerance
                   if (l1(1) < (1-tol) && l1(1) >= -tol && l1(2) < (1-tol) && l1(2) >= -tol) 
                       k2 = k2 + 1;
                       pos2(k2,1) = rr1(1);
                       pos2(k2,2) = rr1(2);
                       pos2(k2,3) = z;
                       id2(k2) = id21;
                       at_num2(k2) = at_num21;
                       lprint2(k2) = true;
                       
                   end
                   if (l2(1) < (1-tol) && l2(1) >= -tol && l2(2) < (1-tol) && l2(2) >= -tol)
                       k2 = k2 + 1;
                       pos2(k2,1) = rr2(1);
                       pos2(k2,2) = rr2(2);
                       pos2(k2,3) = z;
                       id2(k2) = id22;
                       at_num2(k2) = at_num22;
                       lprint2(k2) = true;
                   end
               end
           end
        end

        nat2 = k2;
        nat = nat1 + nat2;
end
