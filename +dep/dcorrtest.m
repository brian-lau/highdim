% DCORRTEST                   Distance correlation test of independence
% 
%     [pval,r,T] = dcorrtest(x,y,varargin)
%
%     Given a sample X1,...,Xn from a p-dimensional multivariate distribution,
%     and a sample Y1,...,Xn from a q-dimensional multivariate distribution,
%     test the hypothesis:
%
%     H0 : X and Y are mutually independent
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL
%     method - 't'          - t-test from Szekely & Rizzo (2013), DEFAULT 
%              'perm'       - randomization using permutation of the rows &
%                             columns of distance matrices
%              'perm-brute' - brute force randomization, directly permuting
%                             one of the inputs, which requires recalculating 
%                             and centering distance matrices
%     nboot - # permutations if not t-test
%
%     OUTPUTS
%     pval - p-value
%     r    - distance correlation, bias-corrected if method = 't' (default)
%     T    - t-statistic
%
%     EXAMPLE
%     rng(1234);
%     p = 2000;
%     n = 500;
%     X = rand(n,p);  Y = X.^2 + 1*randn(n,p);
%
%     tic;[pval,r] = dep.dcorrtest(X,Y,'method','t'); toc
%     [pval , r]
%     tic;[pval,r] = dep.dcorrtest(X,Y,'method','perm','nboot',200);toc
%     [pval , r]
%     tic;[pval,r] = dep.dcorrtest(X,Y,'method','perm-brute','nboot',200);toc
%     [pval , r]
%
%     REFERENCE
%     Szekely et al (2007). Measuring and testing independence by correlation 
%       of distances. Ann Statist 35: 2769-2794
%     Szekely & Rizzo (2013). The distance correlation t-test of independence 
%       in high dimension. J Multiv Analysis 117: 193-213
%
%     SEE ALSO
%     dcorr, DepTest2

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

function [pval,r,T] = dcorrtest(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'method','t',@ischar);
addParamValue(par,'nboot',999,@(x) isnumeric(x) && isscalar(x));
parse(par,x,y,varargin{:});

[n,~] = size(x);
assert(n == size(y,1),'DCORRTEST requires x and y to have the same # of samples');

nboot = par.Results.nboot;

switch lower(par.Results.method)
   case {'t','ttest','t-test'}
      r = dep.dcorr(x,y,'unbiased',true);
      v = n*(n-3)/2;
      T = sqrt(v-1) * r/sqrt(1-r^2);
      pval = 1 - tcdf(T,v-1);
      %tcdf(T,v-1,'upper')
      
      return;
   case {'perm'}      
      a = sqrt(utils.sqdist(x,x));
      b = sqrt(utils.sqdist(y,y));
      [d,dvx,dvy] = dep.dcov(a,b,'dist',true,'unbiased',true);
      r = d/sqrt(dvx*dvy);

      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         [d2,dvx2,dvy2] = dep.dcov(a,b(ind,ind),'dist',true,'unbiased',true);
         boot(i) = d2/sqrt(dvx2*dvy2);
      end
   case {'perm-brute'}
      [d,dvx,dvy] = dep.dcov(x,y,'unbiased',true);
      r = d/sqrt(dvx*dvy);

      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         [d2,dvx2,dvy2] = dep.dcov(x,y(ind,:),'unbiased',true);
         boot(i) = d2/sqrt(dvx2*dvy2);
      end
   otherwise
      error('Unrecognized test method');
end

pval = (1 + sum(boot>r)) / (1 + nboot);
