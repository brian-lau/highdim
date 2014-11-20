function h = rppdf(theta,p)

h = (1/sqrt(pi)) * (gamma(p/2)/(gamma((p-1)/2)))*...
   (sin(theta).^(p-2));
