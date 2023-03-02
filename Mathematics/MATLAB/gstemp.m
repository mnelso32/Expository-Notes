function [Q,R] = gstemp(A)     
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n);
for k = 1:n
   for i = 1:k-1
      R(i,k) = Q(:,i)'*A(:,k);
   end
   tmpAk = A(:,k);
   for i = 1:k-1
      tmpAk = tmpAk - R(i,k)*Q(:,i);
   end
   R(k,k) = norm(tmpAk,2);
   Q(:,k) = tmpAk/R(k,k);
end