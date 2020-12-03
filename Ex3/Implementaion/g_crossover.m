function[pop]=g_crossover(mating_pool,pc,n,m)
    tc=[1 1 1 1 2 2 3 3 1 4];
    
    for tt=1:tc(randi([1 10],1,1))
        par_num=randperm(n);
        for j=1:2:n
            p1=par_num(j);
            p2=par_num(j+1);
            cp=randi([1 m],1,2);
            sort(cp);
            off1=mating_pool(p1,:);
            off2=mating_pool(p2,:);
            if rand<pc
                pos1=[];
                for it=cp(1):cp(2)
                    pos1=[pos1;[find(off2==off1(it)) off1(it)]];
                end
                sortrows(pos1);
                
                pos2=[];
                for it=cp(1):cp(2)
                    pos2=[pos2;[find(off1==off2(it)) off2(it)]];
                end
                sortrows(pos2);
                
                for it=cp(1):cp(2)
                    off1(it)=pos1(it-cp(1)+1,2);
                    off2(it)=pos2(it-cp(1)+1,2);
                end
            end
            pop(j,:)=off1;
            pop(j+1,:)=off2;
        end
    end
 return;