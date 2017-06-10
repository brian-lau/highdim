% HSICTEST                    HSIC test of independence
% 
%     [pval,stat,boot] = hsictest(x,y,varargin)
%
%     Given a sample X1,...,Xm from a p-dimensional multivariate distribution,
%     and a sample Y1,...,Xm from a q-dimensional multivariate distribution,
%     test the hypothesis:
%
%     H0 : X and Y are mutually independent
%
%     INPUTS
%     x - [m x p] m samples of dimensionality p
%     y - [m x q] m samples of dimensionality q
%
%     OPTIONAL (name/value pairs)
%     sigmax - gaussian bandwidth, default = median heuristic
%     sigmay - gaussian bandwidth, default = median heuristic
%     biased - boolean indicated biased estimator (default=false)
%     nboot - # of bootstrap samples, default = 1000
%
%     OUTPUTS
%     pval - p-value
%     stat - HSIC
%     boot - bootstrap samples 
%
%     REFERENCE
%     Gretton et al (2008). A kernel statistical test of independence. NIPS
%     Song et al (2012). Feature Selection via Dependence Maximization.
%       Journal of Machine Learning Research 13: 1393-1434
%
%     SEE ALSO
%     hsic, DepTest2

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

function [pval,stat,boot] = hsictest(x,y,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'x',@isnumeric);
addRequired(par,'y',@isnumeric);
addParamValue(par,'nboot',1000,@(x) isscalar(x)&&isnumeric(x));
parse(par,x,y,varargin{:});

[m,p] = size(x);
[n,q] = size(y);
if m ~= n
   error('x and y must have same # of rows');
end

[stat,K,L,sigmax,sigmay,biased] = dep.hsic(x,y,par.Unmatched);

nboot = par.Results.nboot;
boot = zeros(nboot,1);
for i = 1:nboot
   ind = randperm(m);
   Lb = L(ind,ind);
   boot(i) = dep.hsic_(K,Lb,m,biased);
end

pval = sum(boot>=stat)./nboot;
