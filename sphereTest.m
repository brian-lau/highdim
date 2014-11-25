% SPHERETEST                  Test if data is spherically distributed
% 
%     [pval,stat] = sphereTest(x,varargin)
%
%     Tests whether the covariance matrix of a sample X1, ..., Xn from a 
%     p-dimensional multivariate distribution is proportional to the identity.
%
%     INPUTS
%     x    - [n x p] matrix, n samples with dimensionality p
%
%     OPTIONAL (name/value pairs)
%     test - 'john' - John, Sugiura, Nagao test (JSN)
%            'nagao' - JSN with Box-Bartlett correction
%            'wang' - JSN with correction for large p
%            'sign' - multivariate sign, non-parametric
%            'bcs' - multivariate sign, correction for large p (DEFAULT)
%     approx - use approximation in BCS (default = true)
%
%     OUTPUTS
%     pval - p-value
%     stat - statistic
%
%     REFERENCE
%     Wang, Q and Yao J (2013). On the sphericity test with large-dimensional
%       observations. Electronic Journal of Statistics 7: 2164-2192
%     Zou et al (2014). Multivariate sign-based high-dimensional tests for
%       sphericity. Biometrika 101: 229-236
%
%     SEE ALSO
%     jsn, signtest

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

switch lower(par.Results.test)
   case {'john','nagao','wang'}
      [pval,stat] = jsn(x,par.Results.test);
   case {'sign','s'}
      [pval,stat] = signtest(x,'sign');
   case {'bcs','b'}
      [pval,stat] = signtest(x,'bcs',par.Results.approx);
   otherwise
      error('Unknown test.');
end
