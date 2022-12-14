fun = @efun3;
A = []; b = []; Aeq = []; beq = [];
lb = [-5]; 
ub = [5];
x0 = 1;
nonlcon = @econ3;

x_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global e
for e = 0:0.4:4
    
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    
    f1 = x^2 ;
    f2 = (x-2)^2 ;
  
    x_store = [x_store; x];
    f1_store = [f1_store; f1];
    f2_store = [f2_store; f2];
    exitflag_store = [exitflag_store; exitflag];
    
end

plot(f1_store,f2_store,'*--')

 
e = [0:0.4:4]'; 
VarNames = {'e', 'exit flag', 'x' , 'f1', 'f2'};
T = table(e, exitflag_store, x_store , f1_store ,f2_store , 'VariableNames',VarNames)
