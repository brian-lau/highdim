% SPHERETEST                  Test if data is spherically distributed
% 
%     [pval,stat] = sphereTest(x,varargin)
%
%     INPUTS
%     x    - [n x p] matrix, n samples with dimensionality p
%
%     OPTIONAL
%     test - 
%
%     OUTPUTS
%     pval - p-value
%     stat - corresponding statistic
%
%     REFERENCE
%
%     SEE ALSO
%     uniSphereTest

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     The full license and most recent version of the code can be found on GitHub:
%     https://github.com/brian-lau/spheretest
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
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [pval,stat] = sphereTest(x,varargin)

import sphere.*

if nargin == 2
   [pval,stat] = sphereTest(x,'test',varargin{1});
   return
elseif rem(nargin,2) == 0
   [pval,stat] = sphereTest(x,'test',varargin{1},varargin{2:end});
   return
end

par = inputParser;
addRequired(par,'x',@isnumeric);
addParamValue(par,'test','bcs',@ischar);
addParamValue(par,'approx',true,@(x) isnumeric(x) || islogical(x));
parse(par,x,varargin{:});

[n,p] = size(x);
theta = spatialMedian(x);
U = spatialSign(bsxfun(@minus,x,theta));

% TODO, block process for n large
UtU = U*U';
UtU(sub2ind([n n],1:n,1:n)) = 0;
UtU = sum(UtU(:).^2);

switch par.Results.test
   case {'sign','s'}
      Q = p/n + (n*(n-1)/n^2) * (p/(n*(n-1))) * UtU - 1;
      stat = n*(p+2)*Q/2;
      pval = 1 - chi2cdf(stat,(p+2)*(p-1)/2);
   case {'bcs','b'}
      % Bias-corrected sign test, p = O(n^2)
      Q = (p/(n*(n-1))) * UtU - 1;
      sigma0 = sqrt( 4*(p-1)/(n*(n-1)*(p+2)) );
      
      if par.Results.approx
         % Approximation when x is multivariate normal
         deltanp = n^(-2) + 2*n^(-3);
      else
         % FIXME: not quite working...
         R = sqrt(sum(U.^2,2));
         Rstar = R + U*theta' - sum(theta.^2)./(2*R);
         erk2 = erk(Rstar,2,n);
         deltanp = (1/n^2) * (2 - 2*erk2 + erk2^2) ...
            + (1/n^3) * (8*erk2 - 6*erk2^2 ...
            + 2*erk2*erk(Rstar,3,n) - 2*erk(Rstar,3,n));
      end
      
      stat = (Q - p*deltanp) / sigma0;
      pval = 1 - normcdf(stat);
   otherwise
      error('Unknown test.');
end

function y = erk(Rstar,k,n)
d = sum(1./Rstar)^k;
y = n^(k-1) * sum( Rstar.^(-k)) ./ d;

