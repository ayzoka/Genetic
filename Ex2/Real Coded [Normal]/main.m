function [res,pos]=main(cnt)


n=200;
iter=2000;%hamgra dar nesf, binary dar 3/4

scale=1000;
p_crossover=0.49;
p_mutation=0.1;

m=3;
lb=[0 14.7 0];
ub=[15.1 94.2 5371];

for itr=1:m
    Population(:,itr)=lb(itr)+(ub(itr)-lb(itr))*rand(n,1);
end

[~,F_average(1),F_opt,opt_sol]=F_eval(Population,n,m);
opt_res(1)=F_opt;
itr_res(1)=1;
pos_res=opt_sol;

for itr=2:iter
    p_crossover=min(p_crossover+0.001,0.99);
    p_crossover=max(p_mutation-0.005,0.001);
    [P_select,F_average(itr),F_opt,opt_sol]=F_eval(Population,n,m);
    
    if(F_opt<opt_res(itr-1))
        opt_res(itr)=F_opt;
        itr_res(itr)=itr;
        pos_res=opt_sol;
    else
        opt_res(itr)=opt_res(itr-1);
        itr_res(itr)=itr_res(itr-1);
    end
    
    [mating_pool]=g_roulette_wheel(Population,n,P_select);
    [Population]=g_crossover(mating_pool,p_crossover,n,m,ub,lb);
    [Population]=g_mutation(Population,p_mutation,n,m,scale,ub,lb);

end

res=opt_res(iter);
disp("result: "+num2str(res));
disp(pos_res);
disp(itr_res(iter));
pos=itr_res(iter);

x=1:iter;
f=figure('visible','off');
plot(x,opt_res,'k',x,F_average,'.-k');
xlabel('iteration');
ylabel('function')
legend('optimal','average');
saveas(f,"realcoded"+cnt,'svg');

return;
