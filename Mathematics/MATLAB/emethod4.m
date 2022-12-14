fun = @efun4;
A = []; b = []; Aeq = []; beq = [];
lb = [0,0]; 
ub = [2,2];
x0 = [1,0.5];
nonlcon = @econ4;

x1_store = []; x2_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global e
for e = 0:0.1:1
    
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

 
e = [0:0.1:1]'; 
VarNames = {'e', 'exit flag', 'x1' , 'x2', 'f1', 'f2'};
T = table(e, exitflag_store, x1_store, x2_store, f1_store ,f2_store , 'VariableNames',VarNames)
