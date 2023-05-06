function [H,Q] = HW4_HHrdcUH(A)

% Householder similar reduction of a square matrix A into upper Hessenberg
% 
% Input:  A is the input square matrix, real or complex
% Output: H is an upper Hessenberg, similar to A
%         Q is orthogonal/unitary, such that Q'*A*Q = H numerically
% 
% Copyright (c) 2017, F. Xue
%

[m,n] = size(A);
if m ~= n
    error('A must be an mxn matrix where n >= n.');
end
V = zeros(n,n-2);
for k = 1 : n-2
    xk = A(k+1:n,k);
    signxk1 = sign(xk(1));
    if signxk1 == 0,    signxk1 = 1;    end
    vk = signxk1*norm(xk)*eye(n-k,1)+xk;
    vk = vk/norm(vk);   vk = vk/norm(vk);
    A(k+1:n,k) = -signxk1*norm(xk)*eye(n-k,1);
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - 2*vk*(vk'*A(k+1:n,k+1:n));
    A(:,k+1:n) = A(:,k+1:n) - 2*(A(:,k+1:n)*vk)*vk';
    V(k+1:n,k) = vk;
end

H = A;

Q = eye(n,n);
for k = n-2:-1:1
    Q(k+1:n,k+1:n) = Q(k+1:n,k+1:n)-2*V(k+1:n,k)*(V(k+1:n,k)'*Q(k+1:n,k+1:n));
end

end