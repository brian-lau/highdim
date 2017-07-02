% RBF                         Kernel matrix using Gaussian radial basis
%
%     k = rbf(sigma,x,y,varargin)
%
%     INPUTS
%     sigma - scalar, standard deviation of Gaussian kernel
%     x     - [m x p] m samples of dimensionality p
%
%     OPTIONAL
%     y - [n x p] n samples of dimensionality p
%
%     OUTPUTS
%     k - kernel matrix

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

function k = rbf(sigma,x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'sigma',@isnumeric);
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'approx','none',@ischar);
parse(par,sigma,x,y,varargin{:});

switch lower(par.Results.approx)
   case {'none'}
      % % one-liner scales poorly
      % k = exp(-pdist2(x,y).^2/2*sigma^2);
      xq = sum(x.*x,2);
      if isempty(y)
         k = bsxfun(@plus,xq,xq') - 2*x*x';
      else
         yq = sum(y.*y,2);
         k = bsxfun(@plus,xq,yq') - 2*x*y';
      end
      k = exp(-k/(2*sigma^2));
   case {'random' 'rfm'}
      phix = utils.rfm(sigma,x,par.Unmatched);
      if isempty(y)
         k = phix'*phix;
      else
         phiy = utils.rfm(sigma,y,par.Unmatched);
         k = phix'*phiy;
      end
   otherwise
      error('Unrecognized method');
end
