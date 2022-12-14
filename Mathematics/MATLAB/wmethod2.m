fun = @wfun2;
A = []; b = []; Aeq = []; beq = [];
lb = [-20,-20]; 
ub = [20,20];
x0 = [-2.5,8];
nonlcon = @wcon2;

x1_store = []; x2_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global w
for w = 0:0.1:1
    
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    
    f1 = 2 + (x(1) - 2)^2 + (x(2) - 1)^2 ;
    f2 = 9*x(1) - (x(2) - 1)^2 ;
    
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
