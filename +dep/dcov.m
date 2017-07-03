% DCOV                        Distance covariance
% 
%     [d,dvx,dvy,A,B] = dcov(x,y,varargin)
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL (as name/value pairs, order irrelevant)
%     unbiased - true indicates bias-corrected estimate (default=false)
%     index    - scalar in (0,2], exponent on Euclidean distance, default = 1
%     dist     - true indicates x & y are distance matrices (default=false)
%     doublecenter - true indicates x & y are double-centered distance 
%                matrices (default=false)
%
%     OUTPUTS
%     d   - distance covariance between x & y
%     dvx - x sample distance variance
%     dvy - y sample distance variance
%     A   - double-centered or U-centered distance matrix for x
%     B   - double-centered or U-centered distance matrix for y
%
%     REFERENCE
%     Szekely et al (2007). Measuring and testing independence by correlation 
%       of distances. Ann Statist 35: 2769-2794
%     Szekely & Rizzo (2013). The distance correlation t-test of independence 
%       in high dimension. J Multiv Analysis 117: 193-213
%
%     SEE ALSO
%     dcovtest, dcorr, dcorrtest, rpdcov, fdcov

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

function [d,dvx,dvy,A,B] = dcov(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'unbiased',false,@isscalar);
addParamValue(par,'index',1,@(x) isscalar(x) && (x>0) && (x<=2));
addParamValue(par,'dist',false,@isscalar);
addParamValue(par,'doublecenter',false,@isscalar);
parse(par,x,y,varargin{:});

[n,~] = size(x);
assert(n == size(y,1),'DCOV requires x and y to have the same # of samples');

if par.Results.doublecenter
   % Inputs are already double-centered distance matrices
   A = x;
   B = y;
else
   if par.Results.dist
      % Inputs are euclidean distance matrices
      a = x;
      b = y;
   else
      % Distance matrices, equivalent to pdist2(x,x) and squareform(pdist(x))
      a = sqrt(utils.sqdist(x,x));
      b = sqrt(utils.sqdist(y,y));
   end
   
   if par.Results.index ~= 1
      a = a.^par.Results.index;
      b = b.^par.Results.index;
   end
end

if par.Results.unbiased
   A = utils.ucenter(a);
   B = utils.ucenter(b);
   
   d = sum(sum(A.*B))/(n*(n-3));
   if nargout > 1
      dvx = sum(sum(A.*A))/(n*(n-3));
      dvy = sum(sum(B.*B))/(n*(n-3));
   end
else
   A = utils.dcenter(a);
   B = utils.dcenter(b);
   
   d = sqrt(sum(sum(A.*B))/n^2);
   if nargout > 1
      dvx = sqrt(sum(sum(A.*A))/n^2);
      dvy = sqrt(sum(sum(B.*B))/n^2);
   end
end