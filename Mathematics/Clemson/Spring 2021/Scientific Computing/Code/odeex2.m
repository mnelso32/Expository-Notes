%ode ex2
func = @(t,y) 2.*y;

t=linspace(1,2,100);
yref = 0.5*exp(2.*t);

plot(t,yref,'k--');
hold on;

t=linspace(0,2,1+40);
y = forwardEuler(func, t, 0.5);

plot(t,y,'r-x');

y = backwardEuler(func, t, 0.5);
plot(t,y,'b-o');

legend('exact','forward Euler', 'backward Euler')

hold off;
