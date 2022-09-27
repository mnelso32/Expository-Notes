%first clear the current figure
clf; 

%define the RHS of the ODE: y'=func(t,y)
%Multiplication between terms of y and t need to be done with ".*", same 
% with division. This is because we will be evaluating more than one point
%at a time, so y and t will be vectors and we don't want a vector product
%between t and y when writing t*y for example. .* gives component-wise
%multiplication
func = @(t,y) 2.*y;

%range of t and y to be plotted
trange = [0 3];
yrange = [-2 100];

%this is how many points we are going to plot in t and y direction
nt = 20;
ny = 20;

%calculate the distance between two points in t and y direction
dt = (trange(2)-trange(1)) / nt;
dy = (yrange(2)-yrange(1)) / ny;

%generate the coordinates of all the points on the grid
[T,Y] = meshgrid(trange(1):dt:trange(2), yrange(1):dy:yrange(2));

%evaluate our function at each of those points (giving us the slope)
Z = func(T,Y);

%we need to scale the vector to make them all the same size
L = sqrt(1/dt^2 + Z.^2.*(1/dy^2));
%L = ones(size(z)); %unscaled version

%plot the scaled vector (1, Z) at the given mesh points
%(scaled to fill half of the distance, k.- means blacK, no arrow, solid line)
quiver(T, Y, 1./L, Z./L, 0.5, 'k.-');

% plot more into the same figure
hold on;

%plot some solutions
ezplot('0.5*exp(2.*t)');

%calculate a numerical solution
%a list of points with coordinates (ft, fy) is returned
%[ft fy] = ode45( func, [0 6], 0.5);

%plot the points with crosses, connected with lines, in red
%plot(ft, fy, 'x-r');

% use forward Euler method to compute a solution:
%t = linspace(0,6,20);
%f = forwardEuler( func, t, 0.5);
%plot(t,f,'x-b')


%limit the axes to our range defined above
xlim(trange);
ylim(yrange);

%add labels
xlabel('t');
ylabel('y');

hold off;
