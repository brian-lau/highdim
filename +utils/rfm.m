% RFM                         Random feature maps for Gaussian kernel
%
%     [phi,W] = rfm(sigma,x,varargin)
%
%     INPUTS
%     sigma - scalar, standard deviation of Gaussian kernel
%     x     - [n x d] n samples of dimensionality d
%
%     OPTIONAL
%     sampling - string indicating method for generating random features
%                'uniform' - (DEFAULT)
%                'qmc' - Quasi-Monte Carlo
%                'orf' - Orthogonal Random Features
%                'mm'  - Moment-Matched
%     D        - scalar, dimensionality of feature map
%     W        - [D x d] pre-computed feature map, convenience for a
%                applying feature map to new data
%     complex  - boolean, true returns map as complex
%     skip     - scalar, # initial points to omit (for sampling = 'qmc') 
%     leap     - scalar, # points in between sets (for sampling = 'qmc') 
%     scramble - boolean, scramble sequence (for sampling = 'qmc')
%     state    - scalar, state of qmc generator (for sampling = 'qmc') 
%
%     OUTPUTS
%     phi - feature map
%           [D x n] when 'complex' = true
%           [2D x n] when 'complex' = false, cos and sin components stacked
%     W   - [D x d] feature map
%
%     REFERENCES
%     Felix et al (2016). Orthogonal random features. Advances in Neural 
%       Information Processing Systems, 1975-1983
%     Rahimi & Recht (2007). Random features for large-scale kernel machines.
%       Proc 20th Int Conf on Neural Information Processing Systems, 1177-1184
%     Shen et al (2017). Random features for shift-invariant kernels with 
%       moment matching. Proc 31st AAAI Conf on AI, 2520-2526
%     Sutherland & Schneider (2015). On the error of random fourier features.
%       UAI'15 Proc 31st Conf on Uncertainty in AI, 862-871
%     Yang et al (2014). Quasi-Monte Carlo feature maps for shift-invariant 
%       kernels. Proc 31st Int Conf on Machine Learning (ICML-14), 485-493

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

% TODO
% o ORF should probably be run in blocks
% o SORF

function [phi,W] = rfm(sigma,x,varargin)
persistent pstream; % for qmc

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'sigma',@isnumeric);
addRequired(par,'x',@isnumeric);
addParamValue(par,'sampling','uniform',@ischar);
addParamValue(par,'complex',false,@islogical);
addParamValue(par,'W',[],@ismatrix);
addParamValue(par,'D',2^4,@(x) isnumeric(x) && isscalar(x));
addParamValue(par,'skip',1000,@(x) isnumeric(x) && isscalar(x));
addParamValue(par,'leap',700,@(x) isnumeric(x) && isscalar(x));
addParamValue(par,'scramble',true,@(x) isnumeric(x) || islogical(x));
addParamValue(par,'state',[],@(x) isnumeric(x) && isscalar(x));
parse(par,sigma,x,varargin{:});

[n,d] = size(x);    % # of dimensions
D = par.Results.D;  % # of random bases

if ~isempty(par.Results.W)
   assert(size(par.Results.W,2)==d,'Feature map dimensionality must match input data');
   W = p.Results.W;
else
   switch lower(par.Results.sampling)
      case {'uniform' 'mc'}
         % Random fourier features
         % This version is more accurate, Sutherland & Schneider (2105)
         W = randn(D,d)/sigma;
      case {'mm'}
         if D < d
            warning('Risk of poor approximation for D << d');
         end
         G = randn(D,d);
         W = utils.whiten(G)/sigma;
      case {'qmc'}
         if isempty(pstream) ...
                || ~isa(pstream,'qrandstream') ...
                || (pstream.PointSet.size(2) ~= d)
            pset = haltonset(d,'Skip',par.Results.skip,...
               'Leap',par.Results.leap);
            if par.Results.scramble
               pset = scramble(pset,'RR2');
            end
            % Persistent stream for properly increment draws on subsequent calls
            pstream = qrandstream(pset);
            %fprintf('Halton random stream opened\n')
         end
         
         if ~isempty(par.Results.state)
            pstream.State = par.Results.state;
         end
         %fprintf('Stream state: %g\n',pstream.State);
         omega = pstream.qrand(D);
         W = norminv(omega,0,1)/sigma;
      case {'orf'}
         G = randn(max(d,D),max(d,D));
         [Q,~] = qr(G);
         
         % Chi-distributed with max(d,D) dof
         % X = randn(max(d,D),max(d,D));
         % s = sqrt(sum(X.^2,2));
         s = sqrt(chi2rnd(max(d,D),max(d,D),1));
         % S ensures that the row norms of SQ & G are identically distributed
         S = diag(s);
         
         W = (S*Q)/sigma;
         W = W(1:D,1:d);
      otherwise
         error('Unrecognized sampling method');
   end
end

z = W*x'; % [D x d] * [n x d]'

if par.Results.complex
   phi = (cos(z) - 1i*sin(z)) * sqrt(1/D);
else
   phi = [cos(z) ; sin(z)] * sqrt(1/D);
end
