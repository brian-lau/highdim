% Squared euclidean distance matrix
% Faster than pdist2(x,x) & squareform(pdist(x))
function D = sqdist(X,Y)

assert(all(size(X)==size(Y)),'Input dimensions must match');

if size(X,2) == 1
   Yt = Y';
   XX = X.*X;
   YY = Yt.*Yt;
   D = bsxfun(@plus,XX,YY) - 2*(X*Yt);
   D(D<0) = 0;
   return;
end

Yt = Y';
XX = sum(X.*X,2);
YY = sum(Yt.*Yt);
D = bsxfun(@plus,XX,YY) - 2*(X*Yt);
D(D<0) = 0;