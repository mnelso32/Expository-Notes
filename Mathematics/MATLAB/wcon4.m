function [c,ceq] = wcon4(x) 

c(1) = -x(2) ; 
c(2) = x(2) - x(1); 
c(3) = x(1) + x(2) - 2;
c(4) = (x(1)-1)^2 + (x(2)-0.5)^2-0.25;
c(5) = 0.0125 - (x(1)-1)^2 - (x(2)-0.5)^2;

ceq = []; 
end