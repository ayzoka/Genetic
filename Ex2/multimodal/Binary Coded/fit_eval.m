function [selection_probability,fit,ave_fit,max_fit,opt_sol] = fit_eval(real_val,N,~)
    for(i=1:N)
        x=real_val(i,:);
        f1=1+cos(2*pi*x(1)*x(2));
        f2=floor(exp(-((abs(x(1))+abs(x(2))))/2));
        fit(i)=(f1*f2+x(3));
    end
    selection_probability=1-(fit/sum(fit));
    ave_fit=mean(fit);
    [max_fit,max_loc]=max(fit);
    opt_sol=real_val(max_loc,:);
    return;
end

