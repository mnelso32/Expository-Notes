function y = backwardEuler(func,t,y1)
% compute FE approx to dy/dt = f(t,y) with initial condition y(t1)=y1
% t is a vector of the t points: [t1,t2,...,tn]

% initialize y to be the same size at t
y = 0 * t;

% initial condition
y(1)=y1;

% use BE to find y_i+1 
for i=1:length(y)-1
    % solve y(i+1) - y(i) + dt * f(t(i+1),y(i+1)) = 0
    % define an inline function where y(i+1) is independent variable
    ode_eqn = @(ynext) ynext - y(i) - ( t(i+1) - t(i) ) * func(t(i+1),ynext);
    % solve w/ Matlab's 1d solver, guess y(i) is the solution
    y(i+1) = fzero(ode_eqn,y(i));
end

