function [H,Q,iter] = HW4_QReig(A)

[m,n] = size(A);
if m ~= n
    error('Input matrix must be square!');
end
tol = eps/4;

A = full(A);

tstart1 = tic;
[H,Q] = HW4_HHrdcUH(A);
tend1 = toc(tstart1);
fprintf('Phase I: input matrix reduced to a similar upper Hessenberg in %.2f secs.\n',tend1);

maxiter = n*4;
q = 0;
tmpv = randn(length(A),1);  tmpv = tmpv/norm(tmpv);
if isreal(A) && norm(A*tmpv-(tmpv'*A)')/norm(A,'fro') >= 4*eps
    maxq = n-2;     real_nonsymm = true;
else
    maxq = n-1;     real_nonsymm = false;
end
iter = 1;
tstart2 = tic;
while q < maxq && iter <= maxiter
    
    for k = 1 : n-1
        if abs(H(k+1,k)) <= tol*(abs(H(k,k))+abs(H(k+1,k+1)))
            H(k+1,k) = 0;
        end
    end
    
    oldq = q;
    for j = n-oldq : -1 : 2
        if H(j,j-1) == 0 || (j > 2 && H(j,j-1) ~= 0 && H(j-1,j-2) == 0 && real_nonsymm)
            q = q + 1;
        else
            break;
        end
    end
    subdgH1n2 = diag(H(1:n-q,1:n-q),-1);
    p = find(subdgH1n2 == 0,1,'last');
    if isempty(p),  p = 0;  end
    sizeH22 = n-p-q;
    if q < maxq
        if sizeH22 >= 2
            evs = eig(H(n-q-1:n-q,n-q-1:n-q));
        else
            evs = H(n-q,n-q);
        end
        if isreal(evs) || ~real_nonsymm
            [~,idx] = min(abs(evs-H(n-q,n-q)));
            [H(p+1:n-q,p+1:n-q),GCS] = HW4_SingleShiftedQRstep(H(p+1:n-q,p+1:n-q),evs(idx));
            for k = 1 : sizeH22-1
                Gt = [conj(GCS(1,k)) -GCS(2,k); conj(GCS(2,k)) conj(GCS(1,k))];
                Q(:,p+k:p+k+1) = Q(:,p+k:p+k+1)*Gt;
                H(1:p,p+k:p+k+1) = H(1:p,p+k:p+k+1)*Gt;
                H(p+k:p+k+1,n-q+1:end) = Gt'*H(p+k:p+k+1,n-q+1:end);
            end
        else
            [H(p+1:n-q,p+1:n-q),HVs] = HW4_DoubleShiftedQRstep(H(p+1:n-q,p+1:n-q),H(n-q-1:n-q,n-q-1:n-q));
            for k = 1 : sizeH22-2
                Q(:,p+k:p+k+2) = Q(:,p+k:p+k+2)-Q(:,p+k:p+k+2)*(2*HVs(:,k))*HVs(:,k)';
                H(1:p,p+k:p+k+2) = H(1:p,p+k:p+k+2)-H(1:p,p+k:p+k+2)*(2*HVs(:,k))*HVs(:,k)';
                H(p+k:p+k+2,n-q+1:end) = H(p+k:p+k+2,n-q+1:end)-(2*HVs(:,k))*(HVs(:,k)'*H(p+k:p+k+2,n-q+1:end));
            end
            Q(:,n-q-1:n-q) = Q(:,n-q-1:n-q)-Q(:,n-q-1:n-q)*(2*HVs(1:2,end))*HVs(1:2,end)';
            H(1:p,n-q-1:n-q) = H(1:p,n-q-1:n-q)-H(1:p,n-q-1:n-q)*(2*HVs(1:2,end))*HVs(1:2,end)';
            H(n-q-1:n-q,n-q+1:end) = H(n-q-1:n-q,n-q+1:end)-(2*HVs(1:2,end))*(HVs(1:2,end)'*H(n-q-1:n-q,n-q+1:end));
        end
    end
    fprintf('Iteration %d: %d eigenvalues have converged.\n',iter,q);
    if q >= maxq
        fprintf('All %d eigenvalues have converged in Schur form in %d iterations.\n',n,iter);
        break;
    end
    iter = iter + 1;
end

tend2 = toc(tstart2);

if q < maxq
    fprintf('Maximum iteration (%d) reached; %d eigenvalues have converged.\n',maxiter,q);
end

fprintf('Phase II: the QR iteration took %.2f seconds.\n',tend2);
