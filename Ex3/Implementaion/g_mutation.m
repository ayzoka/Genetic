function [pop]=g_mutation(pop,pm,n,m)

    for it=1:n
        tm=[1 1 1 1 2 2 3 3 4 4];
        for k=1:tm(randi([1 10],1,1))
            if(rand<=pm)
                x=randperm(m,2);
                temp=pop(it,x(1));
                pop(it,x(1))=pop(it,x(2));
                pop(it,x(2))=temp;
            end
        end
    end
    
return;

