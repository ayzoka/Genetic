clc;
tic;
res=[];
t=10;
for it=1:t
    disp("iteration: "+num2str(it))
    clearvars -except res it t

    global gn;
    gn=[1 1 3 4 2 2 3 3 4 1 1];
    N=200;
    Pc=0.99;
    Pm=0.0015;
    ITER=2000; 
    BS=[60 60 60];
    L=sum(BS);
    
    res(it)=genetic(N,Pc,Pm,ITER,BS,L,num2str(it));  
end

runtime=toc;
disp("Crossover array: ");
disp(gn);
mn=min(res);
mx=max(res);
av=sum(res)/t;
disp("average result: "+num2str(av));
disp("minimum result: "+num2str(mn));
disp("maximum result: "+num2str(mx));
disp("running time: "+num2str(runtime));
figure('visible','on')
plot(res);