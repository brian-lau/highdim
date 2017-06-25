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
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL
%     test - 't' indicates t-test from Szekely & Rizzo (2013), 
%              otherwise bootstrap (default = 't')
%     nboot - # bootstrap samples if not t-test
%
%     OUTPUTS
%     pval - p-value
%     r    - distance correlation, corrected if method = 't' (default)
%     T    - t-statistic
%
%     REFERENCE
%     Szekely et al (2007). Measuring and testing independence by correlation 
%       of distances. Ann Statist 35: 2769-2794
%     Huang & Huo (2017). A statistically and numerically efficient
%       independence test based on random projections and distance
%       covariance. arxiv.org/abs/1701.06054v1
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

function [pval,d] = dcovtest(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'estimator','direct',@ischar);
addParamValue(par,'test','perm',@ischar);
addParamValue(par,'nboot',999,@(x) isnumeric(x) && isscalar(x));
parse(par,x,y,varargin{:});

[n,~] = size(x);
if n ~= size(y,1)
   error('DCOVTEST requires x and y to have the same # of samples');
end

switch lower(par.Results.test)
   case {'gamma'}
      % TODO
   otherwise % permutation
      d = dep.dcov(x,y);
      nboot = par.Results.nboot;
      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         boot(i) = dep.dcov(x,y(ind,:));
      end
      pval = (1 + sum(boot>d)) / (1 + nboot);
end