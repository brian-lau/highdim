% Squared euclidean distance matrix
% Faster than pdist2(x,x) & squareform(pdist(x))
function D = sqdist(X,Y)

if (nargin == 1) || isempty(Y)
   XX = sum(X.*X,2);
   D = bsxfun(@plus,XX,XX') - 2*(X*X');
else
   assert(all(size(X)==size(Y)),'Input dimensions must match');
   if size(X,2) == 1
      D = sum((X - Y).^2);
   else
      Yt = Y';
      XX = sum(X.*X,2);
      YY = sum(Yt.*Yt);
      D = bsxfun(@plus,XX,YY) - 2*(X*Yt);
   end
end

D(D<0) = 0;
