function [U,S,V] = HW6_GKLsvds(A,k)

% The Golub-Kahan-Lanczos bidiagonalization
%
% Input:
% A     The matrix for which we are computing largest singular values
%       should be large to show the competitiveness of iterative method
% k     The number of singular values wanted
%
% Output:
% S     The k by k diagonal matrix of approximate dominant singular values
% U,V   The k approximate left and right singular vectors
%
% Copyright (c) F. Xue 10/21/2017

[m,n] = size(A);
v_k = randn(n,1);
v_k = v_k/norm(v_k);
maxiter = min([m n max([ceil(1.2*k) k+5])]);
beta_km1 = 0;
u_km1 = zeros(m,1);
alpha_all = zeros(maxiter,1);
beta_all = zeros(maxiter+1,1);
U = zeros(m,maxiter);
V = zeros(n,maxiter+1);
V(:,1) = v_k;

for iter = 1 : maxiter
    u_k = A*v_k-beta_km1*u_km1;
    
    for jj = 1 : iter-1
        u_k = u_k - U(:,jj)*(U(:,jj)'*u_k);
    end
    
    alpha_k = norm(u_k);
    u_k = u_k/alpha_k;
    v_kp1 = (u_k'*A)'-alpha_k*v_k;
    
    for jj = 1 : iter
        v_kp1 = v_kp1 - V(:,jj)*(V(:,jj)'*v_kp1);
    end
    
    beta_k = norm(v_kp1);
    v_kp1 = v_kp1/beta_k;
    
    alpha_all(iter) = alpha_k;
    beta_all(iter+1) = beta_k;
    
    U(:,iter) = u_k;
    V(:,iter+1) = v_kp1;
    
    u_km1 = u_k;
    v_k = v_kp1;
    beta_km1 = beta_k;
end

B = spdiags([alpha_all beta_all(1:end-1)],0:1,maxiter,maxiter);
[Us,S,Vs] = svd(full(B));
U = U*Us(:,1:k);
V = V(:,1:maxiter)*Vs(:,1:k);
S = S(1:k,1:k);

return
