function F = wfun2(x)
f1 = 2 + (x(1) - 2)^2 + (x(2) - 1)^2 ;
f2 = 9*x(1) - (x(2) - 1)^2 ;
global w;
F = w*f1 + (1-w)*f2 ;
end