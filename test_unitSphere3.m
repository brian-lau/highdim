n = [240];% 1000];
p = 10;%[4 10];%[4 10 20];
kappa = 2;%[0 1 2];%[0 1 2 4];
reps = 1;

tic;
for i = 1:numel(kappa)
   for j = 1:numel(p)
      for k = 1:reps
         x = randvonMisesFisherm(p(j),n,kappa(i))';
         %x = [randn(n,p(j)) ; randvonMisesFisherm(p(j),n,kappa(i))'];
         
%          mu = zeros(1,p(j));
%          mu(end) = 1;
%          x = [randvonMisesFisherm(p(j),n/2,kappa(i),mu)' ;...
%             randvonMisesFisherm(p(j),n/2,kappa(i),-mu)'];
%          mu = zeros(1,p(j));
%          mu(end) = 1;
%          x = [randvonMisesFisherm(p(j),n/3,kappa(i),mu)' ;...
%             randvonMisesFisherm(p(j),n/3,kappa(i),-mu)' ;...
%             randvonMisesFisherm(p(j),n/3,kappa(i),rand(size(mu)))'];

         pval_r(i,j,k) = uniSphereTest(x,'rayleigh');
         pval_rp(i,j,k) = uniSphereTest(x,'rp');
         pval_b(i,j,k) = uniSphereTest(x,'bingham');
         pval_g(i,j,k) = uniSphereTest(x,'gine');
         pval_a(i,j,k) = uniSphereTest(x,'ajne');
         pval_ga(i,j,k) = uniSphereTest(x,'ga');
         [clusteriness, temp, dists, k2] = pairsClusterTest(x);
         pval_p(i,j,k) = temp;
      end
      prob_r(i,j) = sum(squeeze(pval_r(i,j,:))<0.05)/reps;
      prob_rp(i,j) = sum(squeeze(pval_rp(i,j,:))<0.05)/reps;
      prob_b(i,j) = sum(squeeze(pval_b(i,j,:))<0.05)/reps;
      prob_g(i,j) = sum(squeeze(pval_g(i,j,:))<0.05)/reps;
      prob_a(i,j) = sum(squeeze(pval_a(i,j,:))<0.05)/reps;
      prob_ga(i,j) = sum(squeeze(pval_ga(i,j,:))<0.05)/reps;
      prob_p(i,j) = sum(squeeze(pval_p(i,j,:))<0.05)/reps;
   end
   toc
end
