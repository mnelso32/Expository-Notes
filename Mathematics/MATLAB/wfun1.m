function F = wfun1(x)
f1 = x(1);
f2 = x(2);
global w;
F = w*f1+(1-w)*f2;
end