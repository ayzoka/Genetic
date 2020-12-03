clear all;
clc;
n_points=8;
%coordinates of p
%p=[1 6 8;2 8 18;3 4 47;4 39 31;5 25 27;6 38 13;7 33 9;8 26 19];

// p=[[1:n_points]' ,round(rand(n_points,2))];

t_pop=500;

%generate preliminary population
for it=1:t_pop
    x=randperm(n_points);
    j=1;
    while j<it
        if pop(j,:)==x
            j=1;
            x=randperm(n_points);
        else
            j=j+1;
        end    
    end
    pop(it,:)=x;
        
end

%calculate distance of each node
for it=1:t_pop
    sum=0;
    for j=1:(n_points-1)
        sum=sum+sqrt((p(pop(it,j+1),2)-p(pop(it,j),2)).^2);
        sum=sum+((p(pop(it,j+1),3)-p(pop(it,j),3)).^2);
    end
    pop(it,n_points+1)=sum;
end

%sort base on unitage distance (sadly O(n^2) '_')
for it=1:t_pop
    for j=1:it
        if pop(j,n_points+1)>pop(it,n_points+1)
            temp=pop(it,:);
            pop(it,:)=pop(j,:);
            pop(j,:)=temp;
        end
    end
end



pc=0.99;
pm=0.0015;
iter=100*n_points;

tc=round(pc*t_pop);
tm=round(pm*t_pop);
tr=t_pop-(tc+tm);


for it=1:iter

    % selecting population for crossover 
    selected=[];
    for j=1:tr
        selected(j,:)=pop(j,:);
    end  
    
       
    % crossover
    c=randperm(t_pop);
    for j=1:tc
           mating_pool(j,:)=pop(c(j),1:n_points);
    end
  
   
   % crossover of parent1 and parent2
    offspringset=[];
    for z=1:tc
        
        n_parent1=ceil(tc*rand);
        n_parent2=ceil(tc*rand);
        parent1=mating_pool(n_parent1,:);
        parent2=mating_pool(n_parent2,:);
            
        crossoverpoint=ceil((n_points-1)*rand);
        
        parent12=parent1(1:crossoverpoint);
        parent21=parent2(1:crossoverpoint);
        
        s1=[];
        s2=[];
        
        for j=1:numel(parent1)
        
            if sum(parent1(j)~=parent21)==numel(parent21)
                s2=[s2 parent1(j)];
            end
            
            if sum(parent2(j)~=parent12)==numel(parent12)
                s1=[s1 parent2(j)];
            end
            
        end
        
        offspring1=[parent12 s1];
        offspring2=[parent21 s2];
        
        
        % evaluate fitness
        
        offspring1(n_points+1)=0;
        offspring2(n_points+1)=0;
        
        for z=1:(n_points-1)
            offspring1(n_points+1)=offspring1(n_points+1)+sqrt(((p(offspring1(z+1),2)-p(offspring1(z),2)).^2)+(p(offspring1(z+1),3)-p(offspring1(z),3)).^2);
        end
        
        for z=1:(n_points-1)
        offspring2(n_points+1)=offspring2(n_points+1)+sqrt(((p(offspring2(z+1),2)-p(offspring2(z),2)).^2)+(p(offspring2(z+1),3)-p(offspring2(z),3)).^2);
        end        
                
        
        offspringset=[offspringset;offspring1;offspring2];
        
    end
    
    
    % end of crossover
    
    
    % Mutation
    
     % selecting a parent for mutation
    m=randperm(t_pop);
    
    for j=1:tm
        individualsForMutation(j,:)=pop(m(j),1:n_points);
    end
    
   
    %doing the mutation
    Mutatedset=[];
    for z=1:tm
       
        parent=individualsForMutation(z,:);
        
        mutation_pointsoint1=ceil(n_points*rand);
        mutation_pointsoint2=ceil(n_points*rand);
        
        temp=parent(mutation_pointsoint1);
        parent(mutation_pointsoint1)=parent(mutation_pointsoint2);
        parent(mutation_pointsoint2)=temp;
        
        
        % Evaluate fitness of mutateed members
        
        parent(n_points+1)=0;
        
        for z=1:(n_points-1)
            parent(n_points+1)=parent(n_points+1)+sqrt(((p(parent(z+1),2)-p(parent(z),2)).^2)+(p(parent(z+1),3)-p(parent(z),3)).^2);
        end          
        
        Mutatedset=[Mutatedset;parent];
               
        
    end
% End mutation

%new population :

pop2=[selected;offspringset;Mutatedset];

        
   %sort the new population     
    for z=1:size(pop)
        for j=1:z
            if pop(j,n_points+1)>pop(z,n_points+1)
                temp=pop(z,:);
                pop2(z,:)=pop(j,:);
                pop2(j,:)=temp;
            end
        end
    end
    %find the best root
    
    for j=1:t_pop
        pop(j,:)=pop2(j,:);
    end
    
    bestfit(it)=pop(1,n_points+1);
    
    %calculate  mean
    mean=0;
    for j=1:t_pop
       mean=mean+pop(j,n_points+1);
    end
    meanfit(it)=mean/t_pop;    
        
        
end   

best_root=pop(1,1:8)
min_distance=pop(1,9)
%subplot(2,1,1)
%plot(bestfit,'.r')

%subplot(2,1,2)
%plot(meanfit)
    

