% UNISPHERETEST               Test if data is uniformly distributed on unit hypersphere 
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
%     stat - statistic
%
%     REFERENCE
%
%     SEE ALSO
%     sphereTest

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
U = spatialSign(x);

switch lower(par.Results.test)
   case {'rayleigh','r'}
      [pval,stat] = rayleigh(U);
   case {'gine','g'}
      stat = gine(U);
      [pval,boot] = bootstrap('sphere.gine');
   case {'gine3','g3'}
      [pval,stat] = gine3(U);
   case {'ajne','a'}
      stat = ajne(U);
      [pval,boot] = bootstrap('sphere.ajne');
   case {'gine-ajne','ga','ag'}
      stat = gineajne(U);
      [pval,boot] = bootstrap('sphere.gineajne');
   case {'bingham','b'}
      [pval,stat] = bingham(U);
   case {'rp'}
      stat = rp(U,k);
      switch par.Results.dist
         case 'asymp'
            for i = 1:k
               test_cdf = [ stat(:,i) , rpcdf(stat(:,i),p)];
               [~,pval(i)] = kstest(stat(:,i),'CDF',test_cdf);
            end
         otherwise
            Umc = spatialSign(randn(nboot,p));
            u0 = spatialSign(randn(1,p));
            Ymc = acos(Umc*u0');
            for i = 1:k
               [~,pval(i)] = kstest2(stat(:,i),Ymc);
            end
      end
      % TODO add bonferroni correction
      [~,~,adj_p] = fdr_bh(pval,.05,'pdep');
      pval = min(adj_p);
   otherwise
      error('Unknown test.');
end

function [pval,boot] = bootstrap(f)
   boot = zeros(nboot,1);
   for j = 1:nboot
      Umc = sphere.spatialSign(randn(n,p));
      boot(j) = feval(f,Umc);
   end
   pval = sum(boot>=stat)/nboot;
end

end