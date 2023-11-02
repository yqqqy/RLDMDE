clear all; 
clc;
warning('off');
format long;

fun_num = 1;
D=30;
Xmin=-100;
Xmax=100;
ps = 120;
nfe_max=10000*D;

targetbest = [-1400;-1300;-1200;-1100;-1000;-900;-800;-700;-600;-500;-400;-300;
    -200;-100;100;200;300;400;500;600;700;800;900;1000;1100;1200;1300;1400];

fhd=str2func('cec13_func');


[Best_score,Best_pos,cg_curve]=RLDMDE(fhd,ps,D,Xmin,Xmax,nfe_max,fun_num);

disp(Best_score-targetbest(fun_num));

plot(cg_curve)
legend('RLDMDE')
