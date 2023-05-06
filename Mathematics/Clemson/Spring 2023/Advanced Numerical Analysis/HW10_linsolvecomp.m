function HW5_linsolvecomp

load HW10_matrices.mat;

f = randn(length(A),1);
maxiter = 500;
reltol = 1e-8;

tstart1 = tic; 
p = symamd(A); 
PA = speye(size(A));
PA = PA(:,p); 
LA = ichol(A(p,p),struct('type','ict','droptol',1e-3)); 
mfun = @(v) PA*(LA'\(LA\(PA'*v)));
tend1 = toc(tstart1);
fprintf('SPD matrix A: preconditioner setup time %.3f secs.\n',tend1);

tstart2 = tic;
%%% Please replace MATLAB's pcg with your implementation
[x,flagc,~,iterc] = pcg(A,f,reltol,maxiter,mfun);
tend2 = toc(tstart2);
if flagc == 0
    fprintf('P-CG converges in %d iters (%.3f secs) to relative tol %.3e.\n',iterc,tend2,reltol);
else
    fprintf('P-CG fails to converge to relative tol %.3e in %d ters (%.3f secs).\n',reltol,maxiter,tend2);
end

tstart3 = tic;
L = chol(A(p,p),'lower');
x = PA*(L'\(L\(PA'*f)));
tend3 = toc(tstart3);
fprintf('MATLAB''s sparse Cholesky direct solver takes %.3f secs.\n',tend3);
fprintf('nnz of A, incomplete and complete Cholesky factors: %d %d %d.\n',nnz(A),2*nnz(LA),2*nnz(L));


% %% Symmetric indefinite matrix, solved by sqmr (symm. indef preconditioner) or minres (spd preconditioner)
% fprintf('\nSymmetric indefinite B: indefinite preconditioner constructed already.\n');
% f = ones(length(B),1);
% PBT = PB'; LBT = LB';
% mfun = @(v) SB*(PB*(LBT\(DB\(LB\(PBT*(SB*v))))));
% tstart4 = tic;
% [x,flagb,~,iterb] = HW5_sqmr(B,f,reltol,maxiter,mfun); 
% tend4 = toc(tstart4);
% if flagb == 0
%     fprintf('P-SQMR converges in %d iters (%.3f secs) to relative tol %.3e.\n',iterb,tend4,reltol);
% else
%     fprintf('P-SQMR fails to converge to relative tol %.3e in %d ters (%.3f secs).\n',reltol,maxiter,tend4);
% end
% 
% fprintf('\nTransforming the indefinite preconditioner to an SPD preconditioner.\n');
% tstart5 = tic;
% [WB,RB2,~] = HW5_signfctrblkD(DB);
% WBT = WB';
% tend5 = toc(tstart5);
% fprintf('Transformation to SPD preconditioner takes %.3f secs.\n',tend5);
% mfun = @(v) SB*(PB*(LBT\(WB*(RB2\(WBT*(LB\(PBT*(SB*v))))))));
% 
% tstart6 = tic;
% %%% Please replace MATLAB's minres with your implementation
% [x,flagb,~,iterb] = minres(B,f,reltol,maxiter,mfun);
% %[xk,flag,relres,iter,resvec] = HW5_minres(A,b,reltol,maxit,mfun,xk);
% tend6 = toc(tstart6);
% 
% if flagb == 0
%     fprintf('P-MINRES converges in %d iters (%.3f secs) to relative tol %.3e.\n',iterb,tend6,reltol);
% else
%     fprintf('P-MINRES fails to converge to relative tol %.3e in %d ters (%.3f secs).\n',reltol,maxiter,tend6);
% end
% 
% tstart7 = tic;
% [L,D,P,S] = ldl(B);
% x = S*(P*(L'\(D\(L\(P'*(S*f))))));
% tend7 = toc(tstart7);
% fprintf('MATLAB''s sparse LDL direct solver takes %.3f secs.\n',tend7);
% fprintf('nnz of B, incomplete and complete LDL factors: %d %d %d.\n',nnz(B),2*nnz(LB),2*nnz(L));

%% Unsymmetric matrix to be solved by GMRES or bicgstab(ell)
% n = 2000^2;
% A = gallery('neumann', n) + 3*speye(n);
% rng('default');
% dA = sprandn(A);
% A = A+dA;

% find a fill-reducing approximate minimal degree reordering of A
tic; 
p = amd(B); 
telapsed_1 = toc;
P_amd = speye(size(B));
P_amd = P_amd(:,p);
P_amdT = P_amd';

% perform an incomplete sparse LU factorization of A
setup.type = 'ilutp';
setup.droptol = 5e-2;
setup.thresh = 1e-2;
tic; 
[LB,UB,P_ilu] = ilu(B(p,p),setup); 
telapsed_2 = toc;
fprintf('\nUnsymmetric matrix A: preconditioner setup time %.3f secs.\n',telapsed_1+telapsed_2);

% the ILU preconditioner
mfun = @(v)P_amd*(UB\(LB\(P_ilu*(P_amd'*v))));
% the right-hand side vector
f = ones(length(B),1);
% preconditioned bi-conjugate gradient stabilized (2)
tic;
[x,flaga,relres,iter] = bicgstabl(B,f,reltol,ceil(maxiter/4),mfun);
telapsed_3 = toc;
% report the result
if flaga == 0
    fprintf('P-BICGSTABL(2) solve: %d mvps (%.3f secs) to tolerance %.3e.\n',iter*4,telapsed_3,reltol);
else
    fprintf('P-BICGSTABL(2) fails to converge to relative tol %.3e in %d iters (%.3f secs).\n',reltol,maxiter,telapsed_3);
end

% preconditioned induced dimension reduction (4)
ss = 4;
tic;
[X,relres,iter,relresvec] = idr_m(B,f,reltol,ceil(maxiter/4),mfun);
telapsed_4 = toc;
% report the result
if flaga == 0
    fprintf('P-IDR(%d) solve: %d mvps (%.3f secs) to tolerance %.3e.\n',ss,ss+(ss+2)*(iter-1),telapsed_4,reltol);
else
    fprintf('P-IDR(%d) fails to converge to relative tol %.3e in %d iters (%.3f secs).\n',ss,reltol,ss+(ss+2)*(iter-1),telapsed_4);
end

maxit_rst = 60;
maxrst = 10;
tic;
%%% Please replace MATLAB's gmres with your implementation
[x,flaga,relres,iter,resvec] = gmres(@(v)B*mfun(v),f,maxit_rst,reltol,maxrst);
%[x,flaga,relres,iter,resvec] = HW5_gmres(A,f,reltol,maxit_rst,maxrst,mfun,zeros(length(A),1));
telapsed_5 = toc;
if flaga == 0
    fprintf('P-GMRES(%d) solve: %d mvps (%.3f secs) to tolerance %.3e.\n',maxit_rst,(iter(1)-1)*maxit_rst+iter(2),telapsed_5,reltol);
else
    fprintf('P-GMRES(%d) fails to converge to relative tol %.3e in %d iters (%.3f secs).\n',maxit_rst,reltol,maxit_rst*maxrst,telapsed_5);
end

%% finally try exact/complete factorization of A (sparse Gauss elimination)
%fprintf('\nStarting the direct solution of A*x = b.\n');
tic; 
[L,U,P,Q,R] = lu(B); 
telapsed_6 = toc;
% forward and backward substitutions
tic;
x = Q*(U\(L\(P*(R\f))));
telapsed_7 = toc;
% report the result
fprintf('MATLAB''s sparse LU direct solution takes %.3f secs.\n',telapsed_6+telapsed_7);
fprintf('nnz of B, incomplete and complete LU factors: %d %d %d.\n',nnz(B),nnz(LB)+nnz(UB),nnz(L)+nnz(U));

