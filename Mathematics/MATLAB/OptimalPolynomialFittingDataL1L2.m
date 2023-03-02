function [x,p1,l1,p2,l2] = OptimalPolynomialFittingDataL1L2(a,b,deg)

% We assume that length(a)=length(b)>=deg.

m = deg; 
n = length(a);

% Find optimal l1 solution using linprog. The vector p1 returns the coefficients 
% of the polynomial which has optimal l1 distance from b.

beq = b; 
bin = zeros(n,1); 
Aeq = [a.^(m:-1:0), -eye(n)]; 
Ain = [zeros(n,m+1), -eye(n)]; 
c = [zeros(m+1,1); ones(n,1)]; 
[x,l1] = linprog(c,Ain,bin,Aeq,beq); 
p1 = (x(1:m+1));

% Find optimal l2 solution using polyfit. The vector p2 returns the coefficients 
% of the polynomial which has optimal l2 distance from b.

[p2,l2] = polyfit(a,b,m); 