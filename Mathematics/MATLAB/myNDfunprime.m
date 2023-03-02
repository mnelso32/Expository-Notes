function y = myNDfunprime(x)
y(1,1) = 3*x(1)^2;
y(1,2) = -1;
y(2,1) = 1;
y(2,2) = cos(x(2));
end