% RPDCOV                      Randomly projected distance covariance
% 
%     d = rpdcov(x,y,k)
%
%     Estimate (unbiased) distance covariance using Huang & Huo algorithm,
%     which has O(nk (log n + p + q)) complexity and O(n) storage compared
%     to O(n^2(p + q)) complexity and O(n^2) storage of the naive estimator.
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL
%     k - scalar integer, number of random projections, default = 50
%
%     OUTPUTS
%     d - distance covariance between x,y
%     omega_k - distance covariance of k univariate random projections
%
%     REFERENCE
%     Huang & Huo (2017). A statistically and numerically efficient
%        independence test based on random projections and distance
%        covariance. arxiv.org/abs/1701.06054v1
%
%     SEE ALSO
%     fdcov, dcov

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

function [d,omega_k] = rpdcov(x,y,k)

[nx,p] = size(x);
[ny,q] = size(y);

if nx ~= ny
   error('RPDCOV requires x and y to have the same # of samples');
end

% Cp = sqrt(pi) * gamma((p+1)/2) / gamma(p/2);
% Cq = sqrt(pi) * gamma((q+1)/2) / gamma(q/2);
Cp = sqrt(pi) * exp(gammaln((p+1)/2) - gammaln(p/2));
Cq = sqrt(pi) * exp(gammaln((q+1)/2) - gammaln(q/2));

omega_k = zeros(k,1);
for kk = 1:k
   % Project onto random basis on unit hypersphere
   ux = x * sphere.spatialSign(randn(1,p))';
   vy = y * sphere.spatialSign(randn(1,q))';

   % Fast O(n log n) distance covariance
   omega_k(kk) = dep.fdcov(ux,vy);
end

d = mean(Cp*Cq*omega_k);