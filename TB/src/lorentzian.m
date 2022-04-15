function L = lorentzian(x,center,gamma)
   L = 1/pi * (0.5*gamma/((x-center)^2 + (0.5*gamma)^2));
end

