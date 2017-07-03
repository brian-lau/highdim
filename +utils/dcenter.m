% Double-centers distance matrix X
%  X_{ij} - X_{i.}/n - X_{.j}/n + X_{..}/n^2, all i, j
%
% Faster & more memory-efficient than using a centering matrix
% H = eye(n) - ones(n)/n;
% X = H*X*H;
%
function [X,X_j,X__] = dcenter(X)

[n,m] = size(X);
assert(m==n,'DCENTER operates on square, symmetric distance matrices');

X_j = mean(X);
X__ = mean(X_j); % mean(X(:))
X = X - bsxfun(@plus,X_j,X_j') + X__;