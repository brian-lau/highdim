n = 30;
p = [1 2 4 8 16 32 64];
reps = 100;

for i = 1:numel(p)
   for j = 1:reps
      x = rand(30,p(i));
      y = rand(30,p(i));
      r(j,i) = dep.dcorr(x,y);
      
      rstar(j,i) = dep.dcorr(x,y,true);
      T(j,i) = sqrt(n*(n-3)/2-1)*rstar(j,i)/sqrt(1-rstar(j,i)^2);
   end
end

%[pval,dc] =dep.dcorrtest([1 2 3 4 5]',[1.4 1.4 3.5 4.2 4.8]')
% dcor.ttest(c(1,2,3,4,5),c(1.4,1.4,3.5,4.2,4.8))
% 
% 	dcor t-test of independence
% 
% data:  c(1, 2, 3, 4, 5) and c(1.4, 1.4, 3.5, 4.2, 4.8)
% T = 5.6569, df = 4, p-value = 0.002406
% sample estimates:
% Bias corrected dcor 
%            0.942809 

% > dcov(c(1,2,3,4),c(1,1,2,6))
% [1] 1.118034
% > dcor(c(1,2,3,4),c(1,1,2,6))
% [1] 0.8947853
