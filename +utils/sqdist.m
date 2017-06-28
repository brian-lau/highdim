% Squared euclidean distance matrix
% Faster than pdist2(x,x) & squareform(pdist(x))
function D = sqdist(X,Y)

Yt = Y';
XX = sum(X.*X,2);
YY = sum(Yt.*Yt);
D = bsxfun(@plus,XX,YY) - 2*(X*Yt);
D(D<0) = 0;