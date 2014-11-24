
n = 64;
p = [4 8 16 32 48 56 60];

for i = 1:numel(p)
   for j = 1:1000
      x = randn(n,p(i));
%       mu = zeros(1,p(i));
%       sigma = ones(1,p(i));
%       x = mvnrnd(mu,sigma,n);
      %[pm(j,i),statm(j,i)] = sphere.mauchly(x);
      [pj(j,i),statj(j,i)] = sphere.jsn(x);
   end
end


n = [4 8 16 32 64 128 256];
p = [4 8 16 32 64 128 256];
reps = 1000;
% for i = 1:numel(p)
%    for j = 1:numel(n)
%       for k = 1:reps
%          x = randn(n(j),p(i));
%          [pj(i,j,k),statj(i,j,k)] = sphere.jsn(x);
%       end
%    end
%    prob1(i,j) = sum(squeeze(pj(i,j,:))<0.05)/reps;
% end

tic;
for i = 1:numel(n)
   for j = 1:numel(p)
      for k = 1:reps
         x = randn(n(i),p(j));
         pval(i,j,k) = sphere.jsn(x);
      end
      prob5(i,j) = sum(squeeze(pval(i,j,:))<0.05)/reps;
   end
   toc
end
