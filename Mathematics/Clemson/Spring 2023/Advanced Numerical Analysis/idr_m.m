function [xk,relres,iter,resvec] = idr_m(A,b,tol,maxiter,mfun,x0)
% right-sided preconditioned IDR(s) for nonsymmetric linear systems
% Overall a most efficient short recurrence nonsymmetric Krylov solver

% Input parameters
% A             coefficient matrix
% b             right-hand side
% tol           relative tolerance, e.g., 1e-6
% maxiter       maximum number of iterations
% mfun          function handle that implements the action of preconditioning

% Output parameters
% xk            the approximate solution
% relres        relative residual norm ||b-Ax||/||b||
% iter          total number of IDR(s) iterations taken
% resvec        the vector of relative residual norms after each mvp

% Default parameter s = 4 for IDR(s)
s = 4;
n = length(b);
if ismatrix(A) && length(A) > 1
    afun = @(x) A*x;
else
    afun = @(x) A(x);
end
if nargin < 5
    mfun = @(v)v;   % no preconditioner
end
if nargin < 6
    x0 = zeros(n,1);   % no preconditioner
end
rk = b - afun(x0);
rng('default');
[P,~,~] = qr(randn(n,s),0);
xk = zeros(n,1);
dx = zeros(n,s);
dr = zeros(n,s);
resvec = zeros((maxiter+1)*(s+1),1);
resvec(1) = norm(rk);
relres = resvec(1)/norm(b);
mvps = 1;       iter = 1;
for ii = 0 : s-1
    %%%prk = mfun(rk);
    Aprk = afun(mfun(rk));
    mvps = mvps + 1;
    omega = real(Aprk'*rk)/(Aprk'*Aprk);
    dx(:,ii+1) = omega*rk;      xk = xk + dx(:,ii+1);
    dr(:,ii+1) = -omega*Aprk;   rk = rk + dr(:,ii+1);
    resvec(mvps) = norm(rk);
    relres = resvec(mvps)/norm(b);
    if relres <= tol,   break;      end
end
ii = s;         jj = 1;
Ptdr = P'*dr;   Ptrk = P'*rk;
warning('off','MATLAB:nearlySingularMatrix');
while iter <= maxiter && relres > tol
    for kk = 0 : s
        C = Ptdr\Ptrk;
        Aprk = rk - dr*C;   mvps = mvps + 1;
        if kk < 1
            AMApr0 = afun(mfun(Aprk));
            omega = (AMApr0'*Aprk)/(AMApr0'*AMApr0);	
        end
        dx(:,jj) = omega*Aprk - dx*C;
        dr(:,jj) = -afun(mfun(dx(:,jj)));
        rk = rk + dr(:,jj);
        xk = xk + dx(:,jj);
        resvec(mvps) = norm(rk);
        relres = resvec(mvps)/norm(b);      
        if relres <= tol,     break;      end
        dm = P'*dr(:,jj);	Ptdr(:,jj) = dm;    Ptrk = Ptrk + dm;
        ii = ii + 1;        jj = jj + 1;	    jj = mod(jj-1,s)+1;
    end
    iter = 1+(mvps-s)/(s+2);
end
warning('on','MATLAB:nearlySingularMatrix');
xk = x0 + mfun(xk);
rk = b - afun(xk);
resvec = resvec(1:mvps);    relres = norm(rk)/norm(b);
if resvec(end)/norm(b) > tol
    fprintf('IDR(%d) did not converge in %d iterations.\n',s,maxiter);
    fprintf('Final relative res %.2d (tol %.2d).\n',relres,tol);
end
