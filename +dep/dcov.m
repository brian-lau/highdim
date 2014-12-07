% DCOV                        Distance covariance
% 
%     [d,dvx,dvy] = dcov(x,y,modified)
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL
%     correct - boolean indicating bias-correction (default=false)
%
%     OUTPUTS
%     d - distance covariance between x,y
%     dvx - x sample distance variance
%     dvy - y sample distance variance
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

function [d,dvx,dvy] = dcov(x,y,correct)

if nargin < 3
   correct = false;
end

[n,~] = size(x);
if n ~= size(y,1)
   error('DCOV requires x and y to have the same # of samples');
end

a = pdist2(x,x);
b = pdist2(y,y);
A = bsxfun(@minus,a,mean(a));
A = bsxfun(@minus,A,mean(a,2));
A = bsxfun(@plus,A,mean(mean(a)));
B = bsxfun(@minus,b,mean(b));
B = bsxfun(@minus,B,mean(b,2));
B = bsxfun(@plus,B,mean(mean(b)));

if correct
   Astar = (n/(n-1)) * (A - a/n);
   Astar(1:(n+1):n*n) = (n/(n-1)) * (mean(a) - mean(mean(a)));
   Bstar = (n/(n-1)) * (B - b/n);
   Bstar(1:(n+1):n*n) = (n/(n-1)) * (mean(b) - mean(mean(b)));
   
   d = dvmod(Astar,Bstar,n);
   if nargout > 1
      dvx = dvmod(Astar,Astar,n);
      dvy = dvmod(Bstar,Bstar,n);
   end
else
   d = dv(A.*B,n);
   if nargout > 1
      dvx = dv(A.*A,n);
      dvy = dv(B.*B,n);
   end
end

% Sample distance variances
function z = dv(x,n)
z = sqrt(sum(sum(x))/n^2);

function z = dvmod(x,y,n)
U = x.*y;
U(1:(n+1):n*n) = 0;
U = sum(U(:)) - (2/(n-2))*(diag(x)'*diag(y));
z = U/(n*(n-3));
