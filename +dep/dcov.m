function [d,dvx,dvy] = dcov(x,y,modified)

if nargin < 3
   modified = false;
end

[n,~] = size(x);
if n ~= size(y,1)
   error('DCOV requires x and y to have the same # of samples');
end

a = pdist2(x,x);
b = pdist2(y,y);
A = bsxfun(@minus,a,mean(a));
A = bsxfun(@minus,A,mean(a,2));
A = bsxfun(@plus,A,mean(mean(a)));
B = bsxfun(@minus,b,mean(b));
B = bsxfun(@minus,B,mean(b,2));
B = bsxfun(@plus,B,mean(mean(b)));

if ~modified
   d = dv(A.*B,n);
   if nargout > 1
      dvx = dv(A.*A,n);
      dvy = dv(B.*B,n);
   end
else
   Astar = (n/(n-1)) * (A - a/n);
   Astar(1:(n+1):n*n) = (n/(n-1)) * (mean(a) - mean(mean(a)));
   Bstar = (n/(n-1)) * (B - b/n);
   Bstar(1:(n+1):n*n) = (n/(n-1)) * (mean(b) - mean(mean(b)));
   
   d = dvmod(Astar,Bstar,n);
   if nargout > 1
      dvx = dvmod(Astar,Astar,n);
      dvy = dvmod(Bstar,Bstar,n);
   end
end

% Sample distance variances
function z = dv(x,n)
%z = sqrt(sum(sum(x))/n^2);
z = sum(sum(x))/n^2;

function z = dvmod(x,y,n)
U = x.*y;
U(1:(n+1):n*n) = 0;
U = sum(U(:)) - (2/(n-2))*(diag(x)'*diag(y));
z = U/(n*(n-3));
