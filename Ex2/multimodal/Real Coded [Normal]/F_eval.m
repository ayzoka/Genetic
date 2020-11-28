function [P_select,F_avg,F_opt,opt_sol] = F_eval(val,n,~)
    F=zeros(1,n);
    for i=1:n
        x=val(i,:);
        f1=1+cos(2*pi*x(1)*x(2));
        f2=floor(exp(-((abs(x(1))+abs(x(2))))/2));
        F(i)=(f1*f2+x(3));
    end
    P_select=1-(F/sum(F));
    F_avg=mean(F);
    [F_opt,pos]=max(F);
    opt_sol=val(pos,:);
    return;
end

