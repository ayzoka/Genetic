function [res,pos]=main(cnt)
n=50;
iter=200;

pc=0.99;
pm=0.0015;

m=10;
var=randi([1 1000],m,2);
% var=[1 1;0 0;3 2;4 1];

pop=zeros(n,m);
for itr=1:n
    pop(itr,:)=randperm(m);
end

[~,ave(1),F_opt,opt_sol]=F_eval(pop,var,n,m);
opt_res(1)=F_opt;
itr_res(1)=1;
pos_res=opt_sol;

for itr=2:iter
    pc=min(pc+0.001,0.99);
    pm=max(pm-0.005,0.001);
    [P_select,ave(itr),F_opt,opt_sol]=F_eval(pop,var,n,m);
    
    if(F_opt<opt_res(itr-1))
        opt_res(itr)=F_opt;
        itr_res(itr)=itr;
        pos_res=opt_sol;
    else
        opt_res(itr)=opt_res(itr-1);
        itr_res(itr)=itr_res(itr-1);
    end
    
    [mating_pool]=g_roulette_wheel(pop,n,P_select);
    [pop]=g_crossover(mating_pool,pc,n,m);
    [pop]=g_mutation(pop,pm,n,m);

end

res=opt_res(iter);
disp("result: "+num2str(res));
disp(pos_res);
disp(itr_res(iter));
pos=itr_res(iter);

x=1:iter;
f=figure('visible','off');
plot(x,opt_res,'k',x,ave,'.-k');
xlabel('iteration');
ylabel('function')
legend('optimal','average');
saveas(f,"TSP"+cnt,'svg');
% scatter(var(:,1),var(:,2))
return;
