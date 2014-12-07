% Modified Gram-Schmidt orthonormalization
%
%http://ocw.mit.edu/courses/mathematics/18-06-linear-algebra-spring-2010/related-resources/MIT18_06S10_gramschmidtmat.pdf
%http://cavern.uark.edu/~arnold/4353/CGSMGS.pdf
function [Q,R] = mgs(A)

[m,n] = size(A);
V = A;
Q = zeros(m,n);
R = zeros(n,n);
for i = 1:n
   R(i,i) = norm(V(:,i));
   Q(:,i) = V(:,i)/R(i,i);
   for j = i+1:n
      R(i,j) = Q(:,i)'*V(:,j);
      V(:,j) = V(:,j)-R(i,j)*Q(:,i);
   end
end