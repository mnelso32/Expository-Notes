function y=methodB(func,t,y1)
% t = [t1,t2,...,tn] has n points and n-1 steps

n = length(t);
y = 0 * t;
y(1)=y1;

for k=1:length(y)-1              
   h = t(k+1) - t(k);         
   fk = func(t(k),y(k)); 	
   y(k+1) = y(k) + h * func(t(k)+ (h/2),y(k)+(h/2)*fk); 
end;