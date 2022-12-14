function [c,ceq] = econ2(x) 

global e;
c(1) = x(1)^2 + x(2)^2 - 225 ; 
c(2) = x(1) - 3*x(2) + 10 ; 
c(3) = 9*x(1) - (x(2)-1)^2 - e;

ceq = []; 
end

