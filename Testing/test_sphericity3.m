% reproduce table 1 from
% Zou et al (2014). Multivariate sign-based high-dimensional tests for
%   sphericity. Biometrika

n = [40 80];%[20 40 60 80];
p = [55 181 642];%[38 55 89 181 331 642];
reps = 2000;
v = [0 0.125 0.250];

tic;
for i = 1:numel(n)
   for j = 1:numel(p)
      for k = 1:numel(v)
         for m = 1:reps
            
            y = randn(n(i),p(j));
            
            vp = round(v(k)*p(j));
            A = [sqrt(2)*ones(vp,1) ; ones(p(j)-vp,1)];
            x = (diag(A)*y')';
            pval(m) = sphereTest(x,'test','bcs','approx',true);
         end
         prob(i,j,k) = mean(pval<=0.05);
      end
      toc
   end
end

100*prob

% reps = 2000 % 24.11.2014
% ans(:,:,1) =
% 
%     4.7000    5.7500    5.8000
%     6.2500    3.9500    4.7000
% 
% ans(:,:,2) =
% 
%    45.1500   47.8500   49.9000
%    87.6500   93.3000   94.1500
% 
% ans(:,:,3) =
% 
%    64.7000   69.6500   72.4500
%    98.8000   99.3500   99.6500


pZ(:,:,1) = [...
4.9 4.9 5.1;...
4.7 5.2 5.1];
pZ(:,:,2) = [...
41 47 49;...
84 91 94];
pZ(:,:,3) = [...
64 68 72;...
99 100 100]
