function y = forwardEuler(func,t,y1)
% solve the ODE y'=f(t,y)
% input: func is the name of the f function, 
% t is a vector of the t points: [t1,t2,...,tn]
% y1 is the initial condition
% output: the vector y 

% initialize y to be the same size at t, and let N be the total number of
% y points we want to find
y = 0 * t;
N = length(y);

% Set initial condition
y(1)=y1;

% use FE to find y_i+1 
for i=1:N-1
    y(i+1) = y(i) + ( t(i+1) - t(i) ) * func(t(i),y(i));
end

