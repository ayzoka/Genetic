function[Population]=g_crossover(mating_pool,p_crossover,n,m,ub,lb)
    par_index=randperm(n);
    for j=1:2:n
        p1=par_index(j);
        p2=par_index(j+1);
        P1=mating_pool(p1,:);
        P2=mating_pool(p2,:);
        if rand<p_crossover
            a=rand(1,m);
            off1=a.*P1+(1-a).*P2;
            off2=a.*P2+(1-a).*P1;
        else
            off1=P1;
            off2=P2;
        end
        Population(j,:)=off1;
        Population(j+1,:)=off2;
    end
    for i=1:n
        for j=1:m
            if Population(i,j)>ub(j)||Population(i,j)<lb(j)||real(Population(i,j))~=Population(i,j)
                Population(i,j)=lb(j)+(ub(j)-lb(j))*rand();
            end
        end
    end
	return;