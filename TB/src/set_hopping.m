%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to set p-p       %%
%% and p-d hoppings          %%
%% according to rotation     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           %%
%% Written by Valerio Vitale %%
%% February 2021             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = set_hopping(type,Vsig,Vpi,l,m,n,theta,layer)

ct = cos(theta(layer+1) - theta(layer));
st = sin(theta(layer+1) - theta(layer));
s2t = sin(2*(theta(layer+1) - theta(layer)));

switch type
    case "pzpz"
       y = (Vsig - Vpi)*n^2 + Vpi;
    case "pxpx" 
       y = ct*((Vsig -Vpi)*l^2 + Vpi) - st*(Vsig - Vpi)*l*m;
    case "pypy"
       y = ct*((Vsig -Vpi)*m^2 + Vpi) + st*(Vsig - Vpi)*l*m;
    case "pxpz"
       y = (Vsig -Vpi)*l*n;
    case "pypz"
       y = (Vsig -Vpi)*m*n;
    case "pxpy"
       y = ct*(Vsig - Vpi)*l*m + st*((Vsig - Vpi)*l^2 + Vpi);
    case "pypx"
       y = ct*(Vsig - Vpi)*l*m - st*((Vsig - Vpi)*m^2 + Vpi);
    case "pzpx"
       y = ct*(Vsig - Vpi)*n*l - st*(Vsig - Vpi)*n*m;
    case "pzpy"
       y = ct*(Vsig - Vpi)*n*m + st*(Vsig - Vpi)*n*l;
    case "dz2pz"
       y = n*(n^2 - 0.5*(l^2 + m^2))*Vsig + sqrt(3)*n*(l^2+m^2)*Vpi;
    case "dxypz"
       y = sqrt(3)*n*l*m*Vsig - 2*n*l*m*Vpi;
    case "dx2y2pz"
       y = 0.5*sqrt(3)*n*(l^2 - m^2)*Vsig - n*(l^2 - m^2)*Vpi;
    case "pzdxy" %dxy dxz
       y = ct*(sqrt(3)*n*l*m*Vsig - 2*n*l*m*Vpi) + st*(sqrt(3)*l*m^2*Vsig + l*(1-2*m^2)*Vpi);
    case "pzdx2y2" %dx2my2 dz^2 dyz
       y = 0.5*(1+ct^2)*(0.5*sqrt(3)*n*(l^2 - m^2)*Vsig - n*(l^2 - m^2)*Vpi) - sqrt(3)*0.5*st^2*(n*(n^2 - 0.5*(l^2 + m^2))*Vsig + sqrt(3)*n*(l^2+m^2)*Vpi) - 0.5*s2t*(sqrt(3)*m*n^2*Vsig - 2*m*n^2*Vpi);
end

end
