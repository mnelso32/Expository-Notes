function [c,ceq] = wcon1(x) 

c(1) = 1 + 0.1*cos(16*atan(x(1)/x(2))) - x(1)^2 - x(2)^2; 
c(2) = (x(1) - 0.5)^2 + (x(2) - 0.5)^2 - 0.5; 

ceq = []; 
end