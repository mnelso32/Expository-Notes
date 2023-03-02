x(1) = 1; x(2) = 2;
f(1) = 11/2; f(2) = 61/11;
for i = 3:20
  x(i) = i;
  f(i) = 111 - (1130 - 3000/f(i-2))/f(i-1);
end;

plot(x,f,'*--')










fun = @wfun1;
A = []; b = []; Aeq = []; beq = [];
lb = [0,0]; 
ub = [pi,pi];
x0 = [0.8,0.8];
nonlcon = @wcon1;

x1_store = []; x2_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global w
for w = 0:0.1:1
    
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    
    f1 = x(1)
    f2 = x(2)
    
    x1_store = [x1_store; x(1)];
    x2_store = [x2_store; x(2)];
    f1_store = [f1_store; f1];
    f2_store = [f2_store; f2];
    exitflag_store = [exitflag_store; exitflag];
    
end

plot(f1_store,f2_store,'*--')

 
w = [0:0.1:1]'; 
VarNames = {'w', 'exit flag', 'x1' , 'x2', 'f1', 'f2'};
T = table(w, exitflag_store, x1_store, x2_store, f1_store ,f2_store , 'VariableNames',VarNames)

