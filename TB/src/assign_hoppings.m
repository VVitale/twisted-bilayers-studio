function hopping = assign_hoppings(a1,a2,delta,neigh_shell,layer,sym_sign,ir,jr,type1,type2,type3,type4,type5,type6)
hopping = zeros(size(delta,1),1);
tol = 0.3;
if (neigh_shell == 1)
    del4=-(2*a1+a2)/3;
    del5=(a1+2*a2)/3;
    del6=(a1-a2)/3;
    for ih = 1 : length(hopping)
        d4=dot(delta(ih,:),del4)/(norm(delta(ih,:))*norm(del4));
        d5=dot(delta(ih,:),del5)/(norm(delta(ih,:))*norm(del5));
        d6=dot(delta(ih,:),del6)/(norm(delta(ih,:))*norm(del6));
        if(abs(d4)-1>tol && abs(d5)-1>tol && abs(d6)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        
        if(sym_sign(ih) < 0)
            % Only hopping of type 4
            if(abs(d4-1)<tol  || abs(d4+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
            elseif(abs(d6-1)<tol  || abs(d6+1)<tol)
                hopping(ih) = -type4(ir,jr(ih),layer);
            else
                hopping(ih) = 0;
            end
        elseif(sym_sign(ih) > 0)
            % type 4 and type 5
            if(abs(d4-1)<tol || abs(d6-1)<tol  || abs(d4+1)<tol  || abs(d6+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
            elseif(abs(d5-1)<tol || abs(d5+1)<tol)
                hopping(ih) = type5(ir,jr(ih),layer);
            else
                hopping(ih) = 0;
            end
        end
    end
elseif(neigh_shell == 2)
    del1 = a1;
    del2 = a1 + a2;
    del3 = a2;
    for ih = 1 : length(hopping)
        d1=dot(delta(ih,:),del1)/(norm(delta(ih,:))*norm(del1));
        d2=dot(delta(ih,:),del2)/(norm(delta(ih,:))*norm(del2));
        d3=dot(delta(ih,:),del3)/(norm(delta(ih,:))*norm(del3));
        if(ir == jr(ih))
            if(abs(d1-1)<tol || abs(d1+1)<tol)
                hopping(ih) = type1(ir,jr(ih),layer);
            elseif(abs(d2-1)<tol || abs(d2+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
            elseif(abs(d3-1)<tol || abs(d3+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
            else
                hopping(ih) = 0;
            end
        elseif(jr(ih) > ir)
            if(sym_sign(ih) < 0)
                %Hopping of type 1 2 and 3
                if(abs(d1+1)<tol)
                    hopping(ih)  = type1(ir,jr(ih),layer);
                elseif(abs(d1-1)<tol)
                    hopping(ih)  = -type1(ir,jr(ih),layer);
                elseif(abs(d2+1)<tol)
                    hopping(ih)  = type2(ir,jr(ih),layer);
                elseif(abs(d3+1)<tol)
                    hopping(ih)  = -type2(ir,jr(ih),layer);
                elseif(abs(d2-1)<tol)
                    hopping(ih)  = -type3(ir,jr(ih),layer);
                elseif(abs(d3-1)<tol)
                    hopping(ih)  = type3(ir,jr(ih),layer);
                else
                    hopping(ih)  = 0;
                end
            elseif(sym_sign(ih) > 0)
                if(abs(d1-1)<tol || abs(d1+1)<tol)
                    hopping(ih) = type1(ir,jr(ih),layer);
                elseif(abs(d2+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
                elseif(abs(d3+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
                elseif(abs(d2-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
                elseif(abs(d3-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
                else
                    hopping(ih) = 0;
                end
            end
        end
    end
elseif(neigh_shell == 3)
    del7=-2*(a1+2*a2)/3;
    del8=2*(2*a1+a2)/3;
    del9=2*(a2-a1)/3;
    for ih = 1 : length(hopping)
        d7=dot(delta(ih,:),del7)/(norm(delta(ih,:))*norm(del7));
        d8=dot(delta(ih,:),del8)/(norm(delta(ih,:))*norm(del8));
        d9=dot(delta(ih,:),del9)/(norm(delta(ih,:))*norm(del9));
        if(abs(d7)-1>tol && abs(d8)-1>tol && abs(d9)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == 9 && jr(ih) == 6)
            if(abs(d7-1) < tol || abs(d8-1) < tol || abs(d9-1) < tol)
                hopping(ih)  = type6(9,6,layer);
            end
        elseif(ir == 11 && jr(ih) == 6)
            if(abs(d7-1) < tol)
                hopping(ih)  = type6(11,6,layer);
            elseif(abs(d8-1) < tol || abs(d9-1) < tol)
                hopping(ih)  = -0.5*type6(11,6,layer);
            end
        elseif(ir == 10 && jr(ih) == 6)
            if(abs(d8-1) < tol)
                hopping(ih)  = -sqrt(3)/2*type6(11,6,layer);
            elseif(abs(d9-1) < tol)
                hopping(ih)  = sqrt(3)/2*type6(11,6,layer);
            end
        elseif(ir == 9 && jr(ih) == 8)
            if(abs(d7-1) < tol)
                hopping(ih)  = type6(9,8,layer);
            elseif(abs(d8-1) < tol || abs(d9-1) < tol)
                hopping(ih)  = -0.5*type6(9,8,layer);
            end
        elseif(ir == 9 && jr(ih) == 7)
            if(abs(d8-1) < tol)
                hopping(ih)  = -sqrt(3)/2*type6(9,8,layer);
            elseif(abs(d9-1) < tol)
                hopping(ih)  = sqrt(3)/2*type6(9,8,layer);
            end
        elseif(ir == 10 && jr(ih) == 7)
            if(abs(d8-1) < tol || abs(d9-1) < tol)
                hopping(ih)  = 3/4*type6(11,8,layer);
            end
        elseif((ir == 11 && jr(ih) == 7) || (ir == 10 && jr(ih) == 8))
            if(abs(d8-1) < tol)
                hopping(ih)  = sqrt(3)/4*type6(11,8,layer);
            elseif(abs(d9-1) < tol)
                hopping(ih)  = -sqrt(3)/4*type6(11,8,layer);
            end
        elseif(ir == 11 && jr(ih) == 8)
            if(abs(d7-1) < tol)
                hopping(ih)  = type6(11,8,layer);
            elseif(abs(d8-1) < tol || abs(d9-1) < tol)
                hopping(ih)  = 1/4*type6(11,8,layer);
            end
        else
            hopping(ih)  = 0;
        end
    end
end
end

