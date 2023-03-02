function [Y,X] = finite_difference_solution(n,func,ua,ub)

% set step size and create vector X = [x0,x1,...,x{n+1}] where x0=0 and x{n+1}=1

h = 1/(n+1);
X = 0:h:1;


% define alpha, beta, and gamma

alpha = (1/h^2) - (1/2*h); 
beta = -2/(h^2); 
gamma = (1/h^2) + (1/2*h); 


% define vector f = [f1,f2,...,fn]

f = zeros(n,1); 
f(1) = func(1/(n+1))-alpha*ua; 
f(n) = func(n/(n+1))-gamma*ub; 
for k=2:n-1    
   f(k)=func(k/(n+1)); 
end;


% define vector Y = [y1,y2,...,yn]

A = diag(alpha*ones(1,n-1),-1) + diag(beta*ones(1,n)) + diag(gamma*ones(1,n-1),1); 
Y = A\f;


% adjoin ua and ub to Y so that Y = [ua,y1,y2,...,yn,ub]

Y = [ua; Y; ub];