% DCOV                        Distance covariance
% 
%     [d,dvx,dvy] = dcov(x,y,correct)
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

% TODO
% x,y can be distance matrices, to simplify permutations

function [d,dvx,dvy,A,B] = dcov(x,y,correct)

if nargin < 3
   correct = false;
end

[n,~] = size(x);
if n ~= size(y,1)
   error('DCOV requires x and y to have the same # of samples');
end

a = sqrt(utils.sqdist(x,x)); % = pdist2(x,x) and squareform(pdist(x))
b = sqrt(utils.sqdist(y,y));

a_j = mean(a);    ai_ = mean(a,2);
b_j = mean(b);    bi_ = mean(b,2);
a__ = (sum(ai_) + sum(a_j)) / (2*n); % mean(a(:))
b__ = (sum(bi_) + sum(b_j)) / (2*n);

A = a - bsxfun(@plus,a_j,ai_) + a__;
B = b - bsxfun(@plus,b_j,bi_) + b__;

% A = a - bsxfun(@plus,mean(a),mean(a,2)) + mean(a(:));
% B = b - bsxfun(@plus,mean(b),mean(b,2)) + mean(b(:));

if correct
   % Astar & Bstar, section 2.4 Szekely & Rizzo
   A = (n/(n-1)) * (A - a/n);
   A(1:(n+1):n*n) = (n/(n-1)) * (a_j - a__);
   B = (n/(n-1)) * (B - b/n);
   B(1:(n+1):n*n) = (n/(n-1)) * (b_j - b__);
   
   d = dvmod(A,B,n);
   if nargout > 1
      dvx = dvmod(A,A,n);
      dvy = dvmod(B,B,n);
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
