clc;
clear;
res=[];
pos=[];
t=30;
tic;
for it=1:t
    clearvars -except it res pos t;
    disp("iteration: "+num2str(it))
    [res(it),pos(it)]=main(num2str(it));
    
end
runtime=toc;
mn=[min(res) min(pos)];
mx=[max(res) max(pos)];
av=[sum(res)/t sum(pos)/t];
disp("average result: "+num2str(av));
disp("minimum result: "+num2str(mn));
disp("maximum result: "+num2str(mx));
disp("running time: "+num2str(runtime));
figure('visible','on')
plot(res);
