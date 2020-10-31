clc;
tic;
a=0.0;
mn=100000;
num=30;
for itr=1:num
    disp("loop: "+num2str(itr))
    clearvars -except a mn num itr gn N Pc Pm ITER BS L Lo Hi

    global gn;
%    gn=[1 1 1 2 2 2 4 3 5 1 1];
    gn=[1 1 1 1 1 1 1 1 1 1 1];
    N=200;
    Pc=0.99;
    Pm=0.0015;
    ITER=2000; 
    BS=[32 32 32];
    L=sum(BS);
    
    temp=genetic(N,Pc,Pm,ITER,BS,L,num2str(itr));
    a=a+temp;
    if(temp<mn)
        mn=temp;
    end
    
end
disp("")
disp("the average result: "+num2str(a/num));
disp("the mnimum result: "+num2str(mn));
disp("gn array: ");
disp(gn);
cc=toc;
disp("running time:"+num2str(cc));