function F = myfun3(x)
f1 = x^2 ;
f2 = (x-2)^2 ;
global w;
F = w*f1+(1-w)*f2;
end