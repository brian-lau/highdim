% U-center distance matrix
%   a_{ij} - a_{i.}/(n-2) - a_{.j}/(n-2) + a_{..}/((n-1)(n-2)), i \neq j
%   and zero diagonal

function [X,X_j,X__] = ucenter(X)

[n,m] = size(X);
assert(m==n,'UCENTER operates on square, symmetric distance matrices');

X_j = sum(X);
X__ = sum(X_j); % sum(X(:))
X = X - bsxfun(@plus,X_j,X_j')/(n-2) + X__/((n-1)*(n-2));
X(1:(n+1):n*n) = 0;
