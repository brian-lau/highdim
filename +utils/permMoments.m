function [mu,sigma2,skew] = permMoments2(A1,A2)

n = size(A1,1);

[T(1),T2(1),S2(1),T3(1),S3(1),U(1),R(1),B(1)] = useful(A1);
[T(2),T2(2),S2(2),T3(2),S3(2),U(2),R(2),B(2)] = useful(A2);

% First three moments
m1 = prod(T)/n + prod(-T)/(n*(n-1));

m2 = prod(S2)/n + ( prod(T.^2-S2) + 2*prod(T2-S2) + 4*prod(-S2) )/(n*(n-1))...
   + ( 4*prod(2*S2-T2) + 2*prod(2*S2-T.^2) ) / (n*(n-1)*(n-2))...
   + prod(2*T2-6*S2+T.^2) / (n*(n-1)*(n-2)*(n-3));

SP1 = prod(S3)/n;
SP2 = ( 4*prod(-S3+U) + 3*prod(T.*S2-S3) + 6*prod(-S3)...
   + 12*prod(-S3+R) + 6*prod(-S3+B) ) / (n*(n-1));
SP3 = ( 3*prod(-T.*S2+2*S3) + prod(T.^3-3*T.*S2+2*S3)...
   + 12*prod(-T.*S2+2*S3-B) + 12*prod(2*S3-R) + 24*prod(2*S3-R-B)...
   + 6*prod(T.*(T2-S2)+2*S3-2*R) + 24*prod(2*S3-U-R)...
   + 8*prod(T3+2*S3-3*R) ) / (n*(n-1)*(n-2));
SP4 = ( 12*prod(T.*S2-6*S3+2*R+2*B) + 6*prod(T.*(-T2+S2)-6*S3+2*U+4*R)...
   + 3*prod(-T.^3+5*T.*S2-6*S3+2*B) + 12*prod(T.*(-T2+2*S2)-6*S3+3*R+2*B)...
   + 8*prod(-6*S3+2*U+3*R) + 24*prod(-T3-6*S3+U+5*R+B) ) / (n*(n-1)*(n-2)*(n-3));
SP5 = ( 3*prod(T.^3+2*T.*(T2-5*S2) + 24*S3-8*R-8*B)...
   + 12*prod(T.*(T2-2*S2) + 2*T3+24*S3-4*U-16*R-4*B) ) / (n*(n-1)*(n-2)*(n-3)*(n-4));
SP6 = prod(-T.^3-6*T.*(T2-3*S2)-8*T3-120*S3+16*U+72*R+24*B)...
   / (n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5));
m3 = SP1 + SP2 + SP3 + SP4 + SP5 + SP6;

mu = m1;
sigma2 = m2 - m1^2;
skew = (m3 - 3*sigma2*m1 - m1^3) / (sigma2^(3/2));

function [T,T2,S2,T3,S3,U,R,B] = useful(A)
T = trace(A);
T2 = sum(sum(A.^2));
S2 = sum(diag(A.^2));
T3 = sum(sum(A^2.*A));
S3 = sum(diag(A).^3);
U = sum(sum(A.^2.*A));
R = diag(A)'*diag(A*A);
B = diag(A)'*A*diag(A);
