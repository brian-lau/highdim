% RVTEST                      Test RV coefficient of dependence
% 
%     [pval,rv,nrv] = rvtest(x,y)
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OUTPUTS
%     pval - p-value from Pearson type III approximation
%     rv   - RV coefficient
%     nrv  - normalized RV coefficient
%
%     REFERENCE
%     Josse et al (2008). Testing the significance of the RV coefficient.
%       Computational Statistics and Data Analysis 53: 82-91
%
%     SEE ALSO
%     rv, dcorr, dcorrtest

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

function [pval,rv,nrv] = rvtest(x,y)

[n,~] = size(x);
if n ~= size(y,1)
   error('RVTEST requires x and y to have the same # of samples');
end

[rv,xx,yy] = dep.rv(x,y);

% mean
bx = trace(xx)^2/trace(xx^2);
by = trace(yy)^2/trace(yy^2);
mu_rv = sqrt(bx*by)/(n-1);

% variance
tx = (n-1)/((n-3)*(n-1-bx)) * ...
     (n*(n+1)*(sum(diag(xx).^2)/trace(xx^2)) - (n-1)*(bx+2));
ty = (n-1)/((n-3)*(n-1-by)) * ...
     (n*(n+1)*(sum(diag(yy).^2)/trace(yy^2)) - (n-1)*(by+2));
var_rv = (2*(n-1-bx)*(n-1-by))/((n+1)*(n-1)^2*(n-2)) *...
     (1 + ((n-3)/(2*n*(n-1)))*tx*ty);

% Standardized RV coefficient
nrv = (rv - mu_rv)/sqrt(var_rv);

% Skewness estimate for Pearson III approximation
[~,~,sk] = utils.permMoments(xx,yy,n);

if sk >= 0
   pval = 1 - gamcdf(nrv-(-2/sk),4/sk^2,sk/2);
else
   pval = gamcdf(sk/abs(sk)*nrv+2/abs(sk),4/sk^2,abs(sk)/2);
end

end