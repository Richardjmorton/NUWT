FUNCTION mydampedsin_tdp, x,p

return,p[0]+p[1]*sin(2.*!pi/(p[2]*(1+P[6]*x))*x-p[3])*exp(-x*p[5])+x*p[4]

END
