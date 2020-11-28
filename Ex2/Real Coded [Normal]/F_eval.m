function [P_select,F_avg,F_opt,opt_sol] = F_eval(val,n,~)
    F=zeros(1,n);
    for(i=1:n)
        x=val(i,:);
        x1=400*(x(1)^0.9);
        x2=22*((-14.7+x(2))^1.2);
        x3=x(3)+1000;
        F(i)=(x1+x2+x3);
    end
    P_select=1-(F/sum(F));
    F_avg=mean(F);
    [F_opt,pos]=min(F);
    opt_sol=val(pos,:);
    return;
end

