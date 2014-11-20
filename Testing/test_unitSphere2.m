
n = [25 50];% 1000];
p = [5];
reps = 100;
b = 0.0;

tic;
for i = 1:numel(n)
   for j = 1:numel(p)
      for k = 1:reps
         mu = b*ones(1,p(j));
         sigma = eye(p(j));
         x = mvnrnd(mu,sigma,n(i));
         
         %pval_r(i,j,k) = uniSphereTest(x,'rayleigh');
         pval_rp0(i,j,k) = uniSphereTest(x,'rp');
         %pval_rp1(i,j,k) = rp(x,5,1);
         %pval_b(i,j,k) = uniSphereTest(x,'bingham');
         %pval_g(i,j,k) = uniSphereTest(x,'gine');
         %pval_g(i,j,k) = uniSphereTest(x,'gine3');
         %pval_g2(i,j,k) = uniSphereTest(x,'gine');
         %[clusteriness, temp, dists, k2] = pairsClusterTest(x);
         %pval_p(i,j,k) = temp;
      end
      %prob_r(i,j) = sum(squeeze(pval_r(i,j,:))<0.05)/reps;
      prob_rp0(i,j) = sum(squeeze(pval_rp0(i,j,:))<0.05)/reps;
      %prob_rp1(i,j) = sum(squeeze(pval_rp1(i,j,:))<0.05)/reps;
      %prob_b(i,j) = sum(squeeze(pval_b(i,j,:))<0.05)/reps;
      %prob_g(i,j) = sum(squeeze(pval_g(i,j,:))<0.05)/reps;
      %prob_g2(i,j) = sum(squeeze(pval_g2(i,j,:))<0.05)/reps;
      %prob_p(i,j) = sum(squeeze(pval_p(i,j,:))<0.05)/reps;
   end
   toc
end
