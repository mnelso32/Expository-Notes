fun = @myfun3;
A = []; b = []; Aeq = []; beq = [];
lb = [-10]; 
ub = [10];
x0 = 1;
nonlcon = @mycon3;

x_store = []; f1_store = []; f2_store = []; exitflag_store = [];

global w
for w = 0:0.1:1
    
    [x,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    
    f1 = x^2 ;
    f2 = (x-2)^2 ;
  
    x_store = [x_store; x];
    f1_store = [f1_store; f1];
    f2_store = [f2_store; f2];
    exitflag_store = [exitflag_store; exitflag];
    
end

plot(f1_store,f2_store,'*--')

 
w = [0:0.1:1]'; 
VarNames = {'w', 'exit flag', 'x' , 'f1', 'f2'};
T = table(w, exitflag_store, x_store , f1_store ,f2_store , 'VariableNames',VarNames)
