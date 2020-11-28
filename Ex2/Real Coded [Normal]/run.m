clc;
clear;
res=[];
pos=[];
t=50;
tic;
for it=1:t
    clearvars -except it res pos t;
    disp("iteration: "+num2str(it))
    [res(it),pos(it)]=main(num2str(it));
    
end
runtime=toc;
mn=[min(res) min(pos)];
mx=[max(res) max(pos)];
avg=[0 0];
for it=1:t
    avg(1)=avg(1)+res(it);
    avg(2)=avg(2)+pos(it);
end
avg=avg/t;
disp("average result: "+num2str(avg));
disp("minimum result: "+num2str(mn));
disp("maximum result: "+num2str(mx));
disp("running time: "+num2str(runtime));
figure("visible","on");
plot(res);
