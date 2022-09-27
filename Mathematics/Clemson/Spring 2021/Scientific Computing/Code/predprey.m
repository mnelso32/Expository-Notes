% Lotka-Volterra equations

f = @(t,y) [0.125*y(1)-0.01*y(1)*y(2); -y(2)+0.01*y(1)*y(2)];

y0=[80;30];
T=100;
[t,y]=ode45(f, [0,T],y0);

%t=linspace(0,T,1000);
%y=forwardEuler2(f,t,y0);


figure(1)
plot(t,y(:,1),'ro-',t,y(:,2),'bo-')
figure(2)
plot(y(:,1),y(:,2),'bx-')

