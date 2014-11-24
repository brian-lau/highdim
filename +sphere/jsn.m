function [pval,T2,h] = jsn(x)

[n,p] = size(x);

S = cov(x,1);
% % Wang & Yao (2013)
% [~,D] = eig(S);
% l = diag(D);
% lbar = mean(l);
% T2 = ((n*p)/2) * (sum((l-lbar).^2)/p) / lbar^2;
% % Nagao (1973) 3.6
% T2 = ((p^2*n)/2) * trace((S./trace(S) - eye(p)./p)^2);
% Ledoit & Wolf (2002)
U = (1/p)*trace((S/((1/p)*trace(S)) - eye(p))^2);
T2 = n*p/2*U;
% 
% John (1972)
% U = (trace(S^2)) / (trace(S))^2;
% T = (p*U-1)/(p-1);
% T2 = (0.5*n*p)*(p-1)*T;

f = 0.5*p*(p+1) - 1;
if 1
   pval = 1 - chi2cdf(T2,f);
%    if T2>chi2inv(0.95,f)
%       h = 1;
%    else
%       h = 0;
%    end
else
   % From Nagao (1973) theorem 5.1
   ap = (1/12)*(p^3+3*p^2-8*p-12-200/p);
   bp = (1/8)*(-2*p^3-5*p^2+7*p+12+420/p);
   cp = (1/4)*(p^3+2*p^2-p-2-216/p);
   dp = (1/24)*(-2*p^3-3*p^2+p+436/p);
   % From Wang & Yao
   % ap = (1/12)*(p^3+3*p^2-12-200/p);
   % bp = (1/8)*(-2*p^3-5*p^2+7*p-12-420/p);
   % cp = (1/4)*(p^3+2*p^2-p-2-216/p);
   % dp = (1/24)*(-2*p^3-3*p^2+p+436/p);
   
%    Pf = 1 - chi2cdf(T2,f);
%    Pf2 = 1 - chi2cdf(T2,f+2);
%    Pf4 = 1 - chi2cdf(T2,f+4);
%    Pf6 = 1 - chi2cdf(T2,f+6);
%    pval = Pf + (1/n)*(ap*Pf6 + bp*Pf4 + cp*Pf2 + dp*Pf);
   Pf = chi2cdf(T2,f);
   Pf2 = chi2cdf(T2,f+2);
   Pf4 = chi2cdf(T2,f+4);
   Pf6 = chi2cdf(T2,f+6);
   P = Pf + (1/n)*(ap*Pf6 + bp*Pf4 + cp*Pf2 + dp*Pf);
   pval = 1 - P;
end