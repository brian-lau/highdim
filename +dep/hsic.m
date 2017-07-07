% HSIC                        Hilbert-Schmidt Independence Criterion
%
%     [stat,K,L,varargout] = hsic(x,y,varargin)
%
%     Estimate the Hilbert-Schmidt Independence Criterion (HSIC) using
%     gaussian kernels.
%
%     INPUTS
%     x - [n x p] n samples of dimensionality p
%     y - [n x q] n samples of dimensionality q
%
%     OPTIONAL (name/value pairs)
%     unbiased - boolean indicated biased estimator (default=false)
%     gram     - true indicates x & y are Gram matrices (default=false)
%     doublecenter - true indicates x & y are double-centered Gram
%                matrices (default=false)
%
%     Additional name/value pairs are passed through to the kernel function.
%
%     OUTPUTS
%     h - Hilbert-Schmidt Independence Criterion
%     K - [n x n] Gram matrix for x
%     L - [n x n] Gram matrix for y
%
%     EXAMPLE
%     rng(1234)
%     n = 1000; p = 50; q = p;
%     x = rand(n,p);
%     y = x.^2;
%     h = dep.hsic(x,y) % default Gaussian kernel with median heuristic
%
%     % Equivalence between distance covariance & HSIC
%     h = dep.hsic(x,y,'kernel','distance');
%     d = dep.dcov(x,y);
%     [4*h d^2]
%
%     REFERENCE
%     Gretton et al (2008). A kernel statistical test of independence. NIPS
%     Sejdinovic et al (2013). Equivalence of distance-based and RKHS-based
%       statistics in hypothesis testing. Annals of Statistics 41: 2263-2291
%     Song et al (2012). Feature Selection via Dependence Maximization.
%       Journal of Machine Learning Research 13: 1393-1434
%
%     SEE ALSO
%     hsictest

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

function [h,K,L,sigmax,sigmay] = hsic(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'kernel','rbf',@ischar);
addParamValue(par,'unbiased',false,@(x) isnumeric(x) || islogical(x));
addParamValue(par,'gram',false,@isscalar);
addParamValue(par,'doublecenter',false,@isscalar);
parse(par,x,y,varargin{:});

[m,p] = size(x);
[n,q] = size(y);

assert(m == n,'HSIC requires x and y to have the same # of samples');
assert(~(par.Results.doublecenter&&par.Results.unbiased),...
   'Cannot compute unbiased HSIC estimate with double-centered Gram matrices.');

if par.Results.doublecenter
   Kc = x;
   Lc = y;
elseif par.Results.gram
   K = x;
   L = y;
else
   switch lower(par.Results.kernel)
      case {'rbf' 'gauss' 'gaussian'}
         [K,sigmax] = utils.rbf(x,[],par.Unmatched);
         [L,sigmay] = utils.rbf(y,[],par.Unmatched);
      case {'distance'}
         % Sejdinovic et al, pg. 2272, example 15
         K = utils.distkern(x,x);
         L = utils.distkern(y,y);
      otherwise
         
   end
end

if par.Results.unbiased % U-statistic
   K = utils.zerodiag(K);
   L = utils.zerodiag(L);
   
   % l = ones(m,1);
   % h = trace(K*L) + (l'*K*l*l'*L*l)/(n-1)/(n-2) - 2*(l'*K*L*l)/(n-2);
   % h = h/(n*(n-3));
   
   % Equivalent, but faster
   Kc = utils.ucenter(K);
   Lc = utils.ucenter(L);
   h = sum(sum(Kc.*Lc))/(n*(n-3));
else                    % V-statistic
   % h = trace(H*K*H*H*L*H)/n^2;
   
   % Equivalent, but faster
   if ~exist('Kc','var')
      Kc = utils.dcenter(K);
      Lc = utils.dcenter(L);
   end
   h = sum(sum(Kc.*Lc))/n^2;
end