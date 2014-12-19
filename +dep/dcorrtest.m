% DCORRTEST                   Distance correlation test of independence
% 
%     [pval,stat] = dcorrtest(x,y,method)
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
%
%     REFERENCE
%     Szekely et al (2007). Measuring and testing independence by correlation 
%       of distances. Ann Statist 35: 2769-2794
%     Szekely & Rizzo (2013). The distance correlation t-test of independence 
%       in high dimension. J Multiv Analysis 117: 193-213
%
%     SEE ALSO
%     dcorr, dcorrtest

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
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
addParamValue(par,'test','t',@ischar);
addParamValue(par,'nboot',1000,@(x) isnumeric(x) && isscalar(x));
parse(par,x,y,varargin{:});

[n,~] = size(x);
if n ~= size(y,1)
   error('DCORRTEST requires x and y to have the same # of samples');
end

switch lower(par.Results.test)
   case {'t','ttest','t-test'}
      r = dep.dcorr(x,y,true);
      v = n*(n-3)/2;
      T = sqrt(v-1) * r/sqrt(1-r^2);
      pval = 1 - tcdf(T,v-1);
   otherwise % bootstrap unmodified 
      r = dep.dcorr(x,y);
      nboot = par.Results.nboot;
      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         boot(i) = dep.dcorr(x,y(ind,:));
      end
      pval = sum(boot>=r)/nboot;
end