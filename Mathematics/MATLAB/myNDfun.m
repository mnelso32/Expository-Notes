function y = myNDfun(x)
y(1,1) = x(1)^3 - x(2);
y(2,1) = x(1) + sin(x(2)) + 3;
end