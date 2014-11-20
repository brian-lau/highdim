% UNISPHERETEST               Test whether data is uniformly distributed on unit hypersphere 
% 
%     [pval,stat] = uniSphereTest(x,varargin)
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
%     

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

% rayleigh - not consistent against alternatives with zero resultant length
%            Mardia & Jupp, pg 209
%            most powerful invariant test against von mises alternative
% bingham  - antipodially symmetric
%            not consistent alterriatives with E[xx'] = (1/p)*Ip
function [pval,stat] = uniSphereTest(x,varargin)

import sphere.*

if nargin == 2
   [pval,stat] = uniSphereTest(x,'test',varargin{1});
   return
elseif rem(nargin,2) == 0
   [pval,stat] = uniSphereTest(x,'test',varargin{1},varargin{2:end});
   return
end

par = inputParser;
addRequired(par,'x',@isnumeric);
addParamValue(par,'test','rayleigh',@ischar);
addParamValue(par,'nboot',1000,@isnumeric);
addParamValue(par,'k',25,@isnumeric);
addParamValue(par,'dist','boot',@ischar);
parse(par,x,varargin{:});

nboot = par.Results.nboot;
k = par.Results.k;
[n,p] = size(x);
U = spheresign(x);

switch lower(par.Results.test)
   case {'rayleigh','r'}
      stat = rayleigh(U);
      pval = 1 - chi2cdf(stat,p);
   case {'gine','g'}
      stat = gine(U);
      [pval,boot] = bootstrap('sphere.gine');
   case {'gine3','g3'}
      if p ~= 3
         error('Only valid for p = 3');
      end
      stat = gine3(U);
      pval = 1 - sumchi2cdf(stat,3);
   case {'ajne','a'}
      stat = ajne(U);
      [pval,boot] = bootstrap('sphere.ajne');
   case {'gine-ajne','ga','ag'}
      stat = gineajne(U);
      [pval,boot] = bootstrap('sphere.gineajne');
   case {'bingham','b'}
      stat = bingham(U);
      pval = 1 - chi2cdf(stat,((p-1)*(p+2))/2);
   case {'rp'}
      stat = rp(U,k);
      
      switch par.Results.dist
         case 'asymp'
            for i = 1:k
               test_cdf = [ stat(:,i) , sphere.rpcdf(stat(:,i),p)];
               [~,pval(i)] = kstest(stat(:,i),'CDF',test_cdf);
            end
         otherwise
            Umc = spheresign(randn(nboot,p));
            u0 = spheresign(randn(1,p));
            Ymc = acos(Umc*u0');
            for i = 1:k
               [~,pval(i)] = kstest2(stat(:,i),Ymc);
            end
      end
      [~,~,adj_p] = fdr_bh(pval,.05,'pdep');
      pval = min(adj_p);
   otherwise
      error('Unknown test.');
end

function [pval,boot] = bootstrap(f)
   boot = zeros(nboot,1);
   for j = 1:nboot
      Umc = spheresign(randn(n,p));
      boot(j) = feval(f,Umc);
   end
   pval = sum(boot>=stat)/nboot;
end

end