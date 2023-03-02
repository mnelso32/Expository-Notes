function [C,a,b,v] = TrackingPortfolio(I,V,X,cash,gamma)

% I(t) = value of S&P 500 at time t 
% V(i,t) = value of stock i at time t. 
% X(i) = number of shares of stock i in the current TP. 
% cash = either new cash to be invested or cash to be taken out. 
% N = number of stocks 
% T = final time. 
% R(t) = yield of S&P at time t.
% r(i,t) = yield of stock i at time t. 
% C = total value of the current TP at time T.
% regressing r(i,:) against R. 
% p has the form p = bx + a.

sz = size(V); 
N = sz(1); 
T = sz(2);
Vf = V(:,T); 
R = zeros(1,T-1);
r = zeros(N,T-1);
a = zeros(N,1); 
b = zeros(N,1); 
v = zeros(N,1); 
p = zeros(N,2);
C = dot(X,Vf) + cash; 

for t = 1:T-1        
    R(1,t) = log(I(t+1)/I(t));
end; 
for i = 1:N     
      for t = 1:T-1         
         r(i,t) = log(V(i,t+1)/V(i,t));     
      end; 
end; 
for i=1:N     
   p(i,:) = polyfit(R,r(i,:),1);  
   a(i) = p(i,2);     
   b(i) = p(i,1);
   v(i) = Vf(i)/(C-gamma*C);
end; 