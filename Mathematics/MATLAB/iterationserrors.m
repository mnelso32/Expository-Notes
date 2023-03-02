function [iterations, errors] = iterationserrors(g,a);
format long;
x = g(a); 
e = x - 2; 
iterations = [x]; 
errors = [e]; 
for i=1:3 	
	x = g(x); 	
	e = x-2; 	
	iterations = [iterations x]; 	
	errors = [errors e]; 
end

for i=1:4 	
	disp([i iterations(i) errors(i)]);
end        
