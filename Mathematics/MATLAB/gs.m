function [Q,R] = gs(A)     
[m,n] = size(A);
Q = zeros(m,n);    
R = zeros(n,n);
for j = 1:n
   for i = 1:j-1 
      R(i,j) = Q(:,i)'*A(:,j);
   end;
   Q(:,j) = A(:,j);
   for i = 1:j-1 
      Q(:,j) = Q(:,j) - R(i,j)*Q(:,i); 
   end;
   Q(:,j) = Q(:,j) / norm(Q(:,j));  
   R(j,j) = Q(:,j)'*A(:,j);
end