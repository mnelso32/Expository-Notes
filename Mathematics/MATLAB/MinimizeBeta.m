function [x,alpha] = MinimizeAlpha(I,V,X,cash,gamma,epsilon,delta,K);

[C,a,b,v] = TrackingPortfolio(I, V, X, cash , gamma);

N = length(v);

c = [zeros(2*N,1);1;0];

Aeq = [v' , zeros(1,N+2) ; zeros(1,N), ones(1,N), 0, 0]; 
beq = [1; K]; 

Ain1 = [(a.*v)', zeros(1,N), -1, 0;
-(a.*v)', zeros(1,N), -1,0;
(b.*v)', zeros(1,N), 0,1;
-(b.*v)', zeros(1,N), 0,-1];
Ain2 = [diag((1-gamma)*v), diag(-delta), zeros(N,2)];
Ain3 = [diag((gamma-1)*v), diag(epsilon), zeros(N,2)];
Ain = [Ain1; Ain2; Ain3];
bin = [0;0;1;-1;zeros(2*N,1)];

intcon = N+1:2*N;
intcon = intcon'
lb = [zeros(2*N+2,1)];
ub = [Inf(N,1); ones(N,1) ; Inf; Inf];

x = intlinprog(c,intcon,Ain,bin,Aeq,beq,lb,ub);
alpha = dot((a.*v),x(1:N));
