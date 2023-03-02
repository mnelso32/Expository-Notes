function y=methodA(func,t,y1)
% t = [t1,t2,...,tn] has n points and n-1 steps

n = length(t);
y = 0 * t;
y(1)=y1;

for k=1:n-1     
   h = t(k+1) - t(k); 
   ode_eqn = @(ynext) ynext - y(k) - h * func(t(k)+h/2,(ynext + y(k))/2);
   y(k+1) = fzero(ode_eqn,y(k)); 
end;