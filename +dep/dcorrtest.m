% DCORRTEST                   Distance correlation test of independence
% 
%     [pval,stat] = dcorrtest(x,y,method)
%
%     INPUTS
%     x
%     y
%
%     OPTIONAL
%     method - 't' indicates high-dimensional t-test, otherwise bootstrap
%
%     OUTPUTS
%     pval -
%     stat - 
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
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [pval,stat] = dcorrtest(x,y,method)

if nargin < 3
   method = 't';
end
nboot = 1000;

[n,p] = size(x);
if n ~= size(y,1)
   error('DCORRTEST requires x and y to have the same # of samples');
end

switch lower(method)
   case {'t','ttest','t-test'}
      rstar = dep.dcorr(x,y,true);
      v = n*(n-3)/2;
      stat = sqrt(v-1) * rstar/sqrt(1-rstar^2);
      pval = 1 - tcdf(stat,v-1);
   otherwise % bootstrap unmodified
      stat = dep.dcorr(x,y);
      boot = zeros(nboot,1);
      for i = 1:nboot
         ind = randperm(n);
         boot(i) = dep.dcorr(x,y(ind,:));
      end
      pval = sum(boot>=stat)/nboot;
end