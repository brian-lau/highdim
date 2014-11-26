function [pval,dc] = dcorrtest(x,y)

nboot = 1000;

[n,p] = size(x);
if n ~= size(y,1)
   error('DCORRTEST requires x and y to have the same # of samples');
end

[d,dvx,dvy] = dep.dcov(x,y);
dc = d/sqrt(dvx*dvy);

denom = sqrt(dvx*dvy);
boot = zeros(nboot,1);
for i = 1:nboot
   ind = randperm(n);
   
   dboot = dep.dcov(x,y(ind,:));
   boot(i) = dboot/denom;
end

pval = sum(boot>=dc)/nboot;
