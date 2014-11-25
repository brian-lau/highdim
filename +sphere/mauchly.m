function [pval,stat] = mauchly(x)

[n,p] = size(x);

% A = (n-1)*cov(x);
% V = det(A)/((1/p)*trace(A))^p; 
% L = V^(n/2);
% 
[~,D] = eig(cov(x));
l = diag(D);
L = ( (prod(l)^(1/p)) / (sum(l)/p) ) ^ (0.5*p*n);

f = 0.5*p*(p+1) - 1;
if 0
   % 
   lL = -2*log(L);
   pval = 1 - chi2cdf(lL,f);
else
   % Box-Bartlett refinement of asymptotic distribution
   rho = 1 - (2*p^2+p+2)/(6*p*n);
   omega = (p+2)*(p-1)*(p-2)*(2*p^3+6*p^2+3*p+2);
   omega = omega / (288*p^2*n^2*rho^2);
   
   lL = -2*rho*log(L);
   Pf = chi2cdf(lL,f);
   Pf4 = chi2cdf(lL,f+4);
   P = Pf + omega*(Pf4 - Pf);
   pval = max(0,1 - P);
end

if nargout == 2
   stat = lL;
end