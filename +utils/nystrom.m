% NYSTROM                     Nystrom approximation of kernel matrix
% 
%     [k,sigma] = rbf(x,y,varargin)
% 
%     INPUTS
%     X     - [n x p] n samples of dimensionality p
% 
%     OPTIONAL
%     k    - scalar, number of columns to sample
%     rsvd - boolean indicating whether to use randomized SVD
%     tol  - scalar tolerance for truncating small singular values
% 
%     Additional name/value pairs are passed through to function for 
%     estimating the kernel when using an approximation method.
% 
%     OUTPUTS
%     phi   - approximate feature mapped data
%     K     - approximate Gram matrix
%
%     REFERENCES
%     Wang (2015). A Practical Guide to Randomized Matrix Computations with 
%       MATLAB Implementations. https://arxiv.org/abs/1505.07570
%
%     SEE ALSO
%     rsvd

function [phi,K] = nystrom(X,varargin)

par = inputParser;
par.KeepUnmatched = true;
addRequired(par,'X',@isnumeric);
addParamValue(par,'k',[],@(x) isnumeric(x) && isscalar(x));
addParamValue(par,'rsvd',false,@islogical);
addParamValue(par,'tol',[],@(x) isnumeric(x) && isscalar(x));
parse(par,X,varargin{:});

n = size(X,1);
if isempty(par.Results.k)
   k = fix(0.5*n);
else
   k = par.Results.k;
end

ind = randperm(n);
ind = ind(1:k);
C = utils.kernel(X,X(ind,:),par.Unmatched); % C = K(:,ind)
W = C(ind,:);

if par.Results.rsvd
   [U,S] = rsvd(W,par.Unmatched);
else
   [U,S] = svd(W);
end
s = diag(S);
if isempty(par.Results.tol)
   tol = max(size(W)) * eps(norm(s,inf)); % from pinv
else
   tol = par.Results.tol;
end
k = sum(s > tol);
s = 1./sqrt(s(1:k));
UW = bsxfun(@times,U(:,1:k),s');

phi = C*UW;

if nargout == 2
   K = phi*phi';
end