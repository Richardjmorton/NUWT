Function mydampedsin, x, p
return, p[0]+p[1]*sin(2.*!pi/p[2]*x-p[3])*exp(-x*p[5])+x*p[4]


end
