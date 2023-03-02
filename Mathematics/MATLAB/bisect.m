function [y,data] = bisect(a,b,func,tol)
% bisect(a,b,func,tol), uses the bisection method to find a root of func
% in (a,b) within tolerance tol

% evaluate the function at the end points
fa = func(a);
fb = func(b);

% Check that fa and fb have opposite signs (since we implicitly assume this)
if fa*fb >= 0    
	error('The function fails to satisfy f(a)*f(b)<0 for inputs a and b')
end

% keep track of progress of algorithm
data = [0 a b fa fb b-a];

% now we can run the algorithm
it_count = 0;
while (b-a) > tol    
	it_count = it_count+1;    
	% Get the midpoint and evaluate it   
	xnew = (a+b)/2;    
	fnew = func(xnew);  

	% determine which interval to keep    
	if fa * fnew <= 0        
		b  = xnew;        
		fb = fnew;    
	else        
		a  = xnew;        
		fa = fnew;    
	end    
	
	datanew = [it_count a b fa fb b-a];    
	data = [data;datanew];
end

% The solution needs to be a number - choose the midpoint of the final interval
y = (a+b)/2;

end