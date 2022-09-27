% stiff ODE example
% 
% stability requires a small enough time step for explicit methods

func = @(t,y) -y.^2 + 10000*(1-exp(-5*t));
y0=100;
T=0.5;

[t,y] = ode45(func,[0 T],y0);
plot(t,y,'-k')
ylim([0,100])
hold on;
 

n_steps=38;
t=linspace(0,T,1+n_steps);

y = forwardEuler(func, t, y0);
max_FE=max(abs(y))
plot(t,y,'r-o')
 
y = backwardEuler(func, t, y0);
max_BE=max(abs(y))
plot(t,y,'b-x')
% 
legend('ode45', 'forward Euler', 'backward euler','location','se')
hold off;
