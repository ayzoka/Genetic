function [P_select,F_avg,F_opt,opt_sol] = F_eval(pop,var,n,m)
    F=zeros(1,n);
    for it=1:n
        d=@(x,y) sqrt(sum((x-y).^2));
        p=pop(it,:);
        x=var(p(1),:);
        y=var(p(m),:);
        F(it)=d(x,y);
        for j=1:m-1
           x=var(p(j),:);
           y=var(p(j+1),:);
           F(it)=F(it)+d(x,y);
        end
    end
    P_select=1-(F/sum(F));
    F_avg=mean(F);
    [F_opt,pos]=min(F);
    opt_sol=pop(pos,:);
    
    return;


