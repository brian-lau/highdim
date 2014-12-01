function M = zerodiag(M)

[m,n] = size(M);

M(1:(m+1):n*m) = 0;
