% define functions for this problem

func = @(t,y) 10*cos(t)-2*y;
funcref = @(t) -3*exp(-2*t)+2*sin(t)+4*cos(t);


% define n time steps as vector t=t(n) where t=[t1,t2,...,tn,tn+1] where t1=0 and tn+1=1 

t = @(n) 0:2^(-n):1;


% define initial value y1=y(0)=1

y1 = 1;


% create vectors of length n+1 for each method containing approximate solutions

yA = @(n) methodA(func,t(n),y1);
yB = @(n) methodB(func,t(n),y1);
yFE = @(n) forwardEuler(func,t(n),y1);


% create vector of length n+1 containing exact solutions

yref = @(n) funcref(t(n));


% calculate errors for each method using norm(-,Inf)

errorA = @(n) norm(yA(n)-yref(n),Inf);
errorB = @(n) norm(yB(n)-yref(n),Inf);
errorFE = @(n) norm(yFE(n)-yref(n),Inf);


% calculate errorrates for each method from n=4 to n=12

errorratesA=[];
errorratesB=[];
errorratesFE=[];
for n=4:12
	errorrateA = errorA(n+1)/errorA(n);
	errorrateB = errorB(n+1)/errorB(n);
	errorrateFE = errorFE(n+1)/errorFE(n);
	errorratesA = [errorratesA errorrateA];
	errorratesB = [errorratesB errorrateB];
	errorratesFE = [errorratesFE errorrateFE];
end;


format longg
for n=4:9
	disp([2^(-n) errorA(n) errorratesA(n)]);
end

0.0625       0.00297363505250337         0.249996916651961
0.03125      0.000743128606277121          0.24999922925091
0.015625      0.000185748409893716         0.249999806743112
0.0078125      4.64371189869972e-05         0.250000008262771
0.00390625      1.16091365649496e-05          0.2499999957156
0.001953125      2.90227519350594e-06         0.249999958380116






for n=4:9
   disp([2^(-n) errorB(n) errorratesB(n)]);
end

0.0625        0.0046803647542335         0.248717212806086
0.03125       0.00112281854007401         0.249359372054863
0.015625       0.00027492791863537          0.24967956766904
0.0078125      6.80259768404134e-05         0.249839852495094
0.00390625       1.6919231358159e-05         0.249919911990386
0.001953125      4.21896890712148e-06           0.2499599771328







for n=4:9
   disp([2^(-n) errorFE(n) errorratesFE(n)]);
end

0.0625         0.120794793141668         0.498877300734096
0.03125         0.059282752942563         0.499439812045252
0.015625           0.0293702988322         0.499720236721317
0.0078125        0.0146182565026653         0.499860143574655
0.00390625       0.00729271634548834          0.49993007806035
0.001953125       0.00364227288089003         0.499965041753024

















