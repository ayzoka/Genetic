function[new_pop]=g_crossover(mating_pool,Pc,N,L)
    global gn;
    r=gn(round((rand(1)*10)+1));
    while r>0
        r=r-1;
        parent_numb=randperm(N);
        for j=1:2:N
            pointer1=parent_numb(j);
            pointer2=parent_numb(j+1);
            cut_point=randi(1,1,L);
            off1=mating_pool(pointer1,:);
            off2=mating_pool(pointer2,:);
            if rand<Pc
                temp=off2;
                off2(cut_point+1:end)=off1(cut_point+1:end);
                off1(cut_point+1:end)=temp(cut_point+1:end);
            end
            new_pop(j,:)=off1;
            new_pop(j+1,:)=off2;
        end
    end
	return;