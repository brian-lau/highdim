function [d,dvx,dvy] = dcov(x,y)

[n,~] = size(x);
if n ~= size(y,1)
   error('DCOV requires x and y to have the same # of samples');
end

A = pdist2(x,x);
B = pdist2(y,y);

Ac = bsxfun(@minus,A,mean(A));
Ac = bsxfun(@minus,Ac,mean(A,2));
Ac = bsxfun(@plus,Ac,mean(mean(A)));
Bc = bsxfun(@minus,B,mean(B));
Bc = bsxfun(@minus,Bc,mean(B,2));
Bc = bsxfun(@plus,Bc,mean(mean(B)));

d = dv(Ac.*Bc,n);
if nargout > 1
   dvx = sum(sum(Ac.*Ac))/n^2;
   dvy = sum(sum(Bc.*Bc))/n^2;
end

% Sample distance variances
function z = dv(x,n)
z = sqrt(sum(sum(x))/n^2);
