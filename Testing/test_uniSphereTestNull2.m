% Compare to Figure 1 of 
% Paindaveine and Verdebout (2014). On high-dimensional sign tests.
%    Submitted to Bernoulli

n = [30 200 1000];
p = 50:50:1000;
reps = 100;

tic;
for i = 1:numel(n)
   for j = 1:numel(p)
      for k = 1:reps
         x = randn(n(i),p(j));
         pval(i,j,k) = uniSphereTest(x);
         %pval(i,j,k) = gine(x);
         %pval(i,j,k) = bingham(x);
      end
      prob1(i,j) = sum(squeeze(pval(i,j,:))<0.01)/reps;
      prob5(i,j) = sum(squeeze(pval(i,j,:))<0.05)/reps;
   end
   toc
end
