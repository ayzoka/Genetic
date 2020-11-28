function [Population]=g_mutation(new_pop,p_mutation,n,m,scale,ub,lb)
	L=ub-lb;
    S=L*scale;
    D=zeros(n,m);
    for j=1:m
        D(:,j)=S(j)*abs(randn(n,1));
    end
    mask= (rand(n,m)<=p_mutation);
    Population=new_pop+(mask.*D);
    for i=1:n
        for j=1:m
            if new_pop(i,j)>ub(j)||new_pop(i,j)<lb(j)||Population(i,j)~=real(Population(i,j))
                new_pop(i,j)=lb(j)+(ub(j)-lb(j))*rand;
            end
        end
    end
return;

