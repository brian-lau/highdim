function [pval,stat] = dcorrtest(x,y,method)

if nargin < 3
   method = 't';
end
nboot = 1000;

[n,p] = size(x);
if n ~= size(y,1)
   error('DCORRTEST requires x and y to have the same # of samples');
end

switch lower(method)
   case {'t','ttest','t-test'}
      rstar = dep.dcorr(x,y,true);
      v = n*(n-3)/2;
      stat = sqrt(v-1) * rstar/sqrt(1-rstar^2);
      pval = 1 - tcdf(stat,v-1);
   otherwise % bootstrap unmodified
      stat = dep.dcorr(x,y);
      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         boot(i) = dep.dcorr(x,y(ind,:));
      end
      pval = sum(boot>=stat)/nboot;
end