func = @(x) -2 - 4*pi*(x+1)*sin(2*pi*x^2) - 16*pi^2*x^2*cos(2*pi*x^2);
funcref = @(x) cos(2*pi.*x.*x)-2.*x;
n=10;
ua=1;
ub=-1;
[Y,X] = finite_difference_solution(n,func,ua,ub);
x = linspace(0,1);
plot(X,Y,x,funcref(x));


x=linspace(0,1);
t=(1/11):(1/11):(10/11);
plot(x,funcref(x),t,u);