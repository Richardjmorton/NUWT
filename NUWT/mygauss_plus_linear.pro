function mygauss_plus_linear,x,p

return, p[0]*exp(-((x-p[1])/p[2])^2/2.)+p[3]+p[4]*x

end
