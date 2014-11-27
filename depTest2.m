% DEPTEST2                    2-sample test of independence
% 
%     [pval,stat] = depTest2(x,y,varargin)
%
%     Tests whether the covariance matrix of a sample X1, ..., Xn from a 
%     p-dimensional multivariate distribution is proportional to the identity.
%
%     INPUTS
%     x    - [n x p] matrix, n samples with dimensionality p
%
%     OPTIONAL (name/value pairs)
%     test - 'dcorr' - Szekely & Rizzo distance correlation 
%            'rv' - RV coefficient
%
%     OUTPUTS
%     pval - p-value
%     stat - statistic
%
%     EXAMPLE
%
%     REFERENCE

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

function [pval,stat] = depTest2(x,y,varargin)

import dep.*

if nargin == 3
   [pval,stat] = depTest2(x,y,'test',varargin{1});
   return
elseif rem(nargin,2) == 1
   [pval,stat] = depTest2(x,y,'test',varargin{1},varargin{2:end});
   return
end

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'test','dcorr',@ischar);
parse(par,x,y,varargin{:});

switch lower(par.Results.test)
   case {'dcorr'}
      [pval,stat] = dcorrtest(x,y);
   case {'rv'}
      [pval,stat] = rvtest(x,y);
   otherwise
      error('Unknown test.');
end
