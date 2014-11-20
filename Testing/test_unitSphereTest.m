% Compare to Figure 1 of Paindaveine and Verdebout

n = [30 200];% 1000];
p = [2 4 8 16 32 64 128 256];% 512 1024];
reps = 100;

tic;
for i = 1:numel(n)
   for j = 1:numel(p)
      for k = 1:reps
         x = randn(n(i),p(j));
         %pval(i,j,k) = rayleigh(x);
         %pval(i,j,k) = gine(x);
         pval(i,j,k) = bingham(x);
      end
      prob1(i,j) = sum(squeeze(pval(i,j,:))<0.01)/reps;
      prob5(i,j) = sum(squeeze(pval(i,j,:))<0.05)/reps;
   end
   toc
end

% Weird failure with symmetric clusters?
mu = [0 0 0];
SIGMA = [1 0 0; 0 1.5 0; 0 0 1];
x = mvnrnd(mu,SIGMA,100);

for i = 1:100
   x = randn(100,5);
   [pval(i),Fn(i)] = gine(x);
end