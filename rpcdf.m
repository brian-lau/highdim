function c = rpcdf(theta,p)

x = 0:.001:pi;
h = rppdf(x,p);

c = cumtrapz(x,h);
c = interp1(x,c,theta);