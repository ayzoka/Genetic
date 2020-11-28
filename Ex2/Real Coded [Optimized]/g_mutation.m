function [res]= g_mutation (parent, mu,LB,UB)
    nVar=numel(parent);
    nmu=ceil(mu*nVar);
    while nmu>0
        nmu=nmu-1;
        j=randsample(nVar,1);
        res=parent;
        res(j)=parent(j)- rand*parent(j);
    end
end