% DCOVTEST                    Distance covariance test of independence
% 
%     [pval,r,T] = dcovtest(x,y,varargin)
%
%     Given a sample X1,...,Xn from a p-dimensional multivariate distribution,
%     and a sample Y1,...,Xn from a q-dimensional multivariate distribution,
%     test the hypothesis:
%
%     H0 : X and Y are mutually independent
%
%     This hypothesis is tested using several different permutation methods. 
%
%     The default permutation method actually avoids permuting the data
%     altogether by approximating the permutation distribution using a 
%     moment-matched Pearson Type III distribution (Bilodeau & Guetsop Nangue
%     2017; Josse et al 2008; Minas & Montana 2014). The first three moments 
%     of the permutation distribution can be calculated exactly for distance 
%     covariance and related statistics (Kazi-Aoual et al 1995), and are 
%     accurate (Josse et al 2008). Since this method does not actually permute 
%     the data, it is very fast, and acheives the same statistical power that 
%     would otherwise require millions of permutations (Minas & Montana, 2014).
%
%     Testing using actual permutations of the data are also implemented.
%     Naive permutation of the rows of X or Y is expensive due to O(n^2) 
%     distance calculations. This can be avoided since it is equivalent to 
%     permuting the rows and columns of the distance matrix simultaneously, 
%     and recomputing the statistic with the permuted distance matrix.
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL (as name/value pairs, order irrelevant)
%     method - 'pearson'    - Pearson type III approx by moment matching (DEFAULT)
%              'perm'       - randomization using permutation of the rows &
%                             columns of the double-centered distance matrices
%              'perm-dist'  - randomization using permutation of the rows &
%                             columns of distance matrices
%              'perm-brute' - brute force randomization, directly permuting
%                             one of the inputs, which requires recalculating 
%                             and centering distance matrices
%     nboot - # bootstrap samples if not t-test
%
%     OUTPUTS
%     pval - p-value
%     d    - distance covariance
%     stat - test statistic
%
%     EXAMPLE
%     rng(1234);
%     p = 100;
%     n = 1000;
%     X = rand(n,p);  Y = X.^2 + 5.5*randn(n,p);
%
%     tic;[pval,d] = dep.dcovtest(X,Y,'method','pearson'); toc
%     pval
%     tic;[pval,d] = dep.dcovtest(X,Y,'method','perm','nboot',200);toc
%     pval
%     tic;[pval,d] = dep.dcovtest(X,Y,'method','perm-brute','nboot',200);toc
%     pval
%
%     REFERENCE
%     Bilodeau & Guetsop Nangue (2017). Approximations to permutation tests 
%       of independence between two random vectors. 
%       Computational Statistics and Data Analysis, submitted.
%     Josse, Pages & Husson (2008). Testing the significance of the RV 
%       coefficient. Computational Statistics and Data Analysis. 53: 82-91
%     Kazi-Aoual et al (1995). Refined approximations to permutation tests 
%       for multivariate inference. Computational Statistics and Data Analysis.
%       20: 643-656
%     Minas & Montana (2014). Distance-based analysis of variance: 
%       Approximate inference. Statistical Analysis & Data Mining. 7: 450-470
%
%     SEE ALSO
%     dcov, dcorr, dcorrtest, DepTest2

%     $ Copyright (C) 2017 Brian Lau, brian.lau@upmc.fr $
%     The full license and most recent version of the code can be found at:
%     https://github.com/brian-lau/highdim
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.

function [pval,d,stat] = dcovtest(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'method','pearson',@ischar);
addParamValue(par,'nboot',999,@(x) isnumeric(x) && isscalar(x));
parse(par,x,y,varargin{:});

[n,~] = size(x);
assert(n == size(y,1),'DCOVTEST requires x and y to have the same # of samples');

nboot = par.Results.nboot;

switch lower(par.Results.method)
   case {'pearson'}
      [d,~,~,A,B] = dep.dcov(x,y);
      [mu,sigma2,skew] = utils.permMoments(A,B); % Exact moments
      
      stat = (n^2)*d^2; %  = sum(sum(A.*B)) for biased estimator
      stat = (stat - mu)/sqrt(sigma2);
      if skew >= 0
         pval = gamcdf(stat - (-2/skew),4/skew^2,skew/2,'upper');
      else
         as = abs(skew);
         pval = gamcdf(skew/as*stat + 2/as,4/skew^2,as/2,'upper');         
      end
   case {'perm'}
      [d,~,~,A,B] = dep.dcov(x,y);

      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         B = B(:,ind); B = B(ind,:);
         boot(i) = dep.dcov(A,B,'doublecenter',true);
      end
      pval = (1 + sum(boot>d)) / (1 + nboot);
   case {'perm-dist'}
      a = sqrt(utils.sqdist(x,x));
      b = sqrt(utils.sqdist(y,y));
      d = dep.dcov(a,b,'dist',true);
      
      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         b = b(:,ind); b = b(ind,:);
         boot(i) = dep.dcov(a,bstar,'dist',true);
      end
      pval = (1 + sum(boot>d)) / (1 + nboot);
   case {'perm-brute'}
      d = dep.dcov(x,y);

      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         boot(i) = dep.dcov(x,y(ind,:));
      end
      pval = (1 + sum(boot>d)) / (1 + nboot);
   otherwise
      error('Unrecognized test method');
end