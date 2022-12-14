fun = @efun2;
A = []; b = []; Aeq = []; beq = [];
lb = [-20,-20]; 
ub = [20,20];
x0 = [-2.5,8];
nonlcon = @econ2;

x1_store = []; x2_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global e
for e = -200:10:-100
    
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

 
e = [-200:10:-100]'; 
VarNames = {'e', 'exit flag', 'x1' , 'x2', 'f1', 'f2'};
T = table(e, exitflag_store, x1_store, x2_store, f1_store ,f2_store , 'VariableNames',VarNames)
