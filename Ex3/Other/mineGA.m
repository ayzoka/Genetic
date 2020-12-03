clear all;
clc;
Ncities=8;
[1  6   8 ;2 8 18 ; 3  4  47 ; 4 39  31 ;5 25 27 ;6  38 13 ; 7  33  9; 8 26  19];
%coordinates of cities

cities=[1  60   81 ;2 84 18 ; 3  40  47 ; 4 39  31 ;5 25 27 ;6  38 130 ; 7  33  90; 8 26  109]
%for i=1:Ncities    
%cities(i,1)=i;
%cities(i,2)=round(80*rand);
%cities(i,3)=round(80*rand);
%end

popsize=500;

citiescoordinates=cities ;


%produce the first population

for i=1:popsize;
    x=randperm(Ncities);
    j=1;
    while j<i;
        if pop(j,:)==x;
            j=1;
            x=randperm(Ncities);
        else
            j=j+1;
        end    
    end
    pop(i,:)=x;
        
end

%pop%first generation



%count the distance of each root
for i=1:popsize;
    x=0;
    for j=1:(Ncities-1);
        x=x+sqrt((cities(pop(i,j+1),2)-cities(pop(i,j),2)).^2)+((cities(pop(i,j+1),3)-cities(pop(i,j),3)).^2);
    end
    pop(i,Ncities+1)=x;
end



%pop%counted the distances

%sort unitage distance

for i=1:popsize;
    for j=1:i
        if pop(j,Ncities+1)>pop(i,Ncities+1);
            temp=pop(i,:);
            pop(i,:)=pop(j,:);
            pop(j,:)=temp;
        end
    end
end



pc=.8;
pm=.08;
iter=100;

numberofcrossover=round(pc*popsize);
numberofmutation=round(pm*popsize);
numberofrec=popsize-(numberofcrossover+numberofmutation);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:iter
    % selecting population for crossover 
    
    selected=[];
    for j=1:numberofrec
        selected(j,:)=pop(j,:);
        
    end  
    
       
    %% crossover
    
    
    c=randperm(popsize);
    for j=1:numberofcrossover;
           individualsForCross(j,:)=pop(c(j),1:Ncities);
    end
  
   
   % crossover of parent1 and parent2
    
    offspringset=[];
    
    for z=1:numberofcrossover;
    
    n_parent1=ceil(numberofcrossover*rand);
    n_parent2=ceil(numberofcrossover*rand);
    parent1=individualsForCross(n_parent1,:);
    parent2=individualsForCross(n_parent2,:);
        
    crossoverpoint=ceil((Ncities-1)*rand);
    
    parent12=parent1(1:crossoverpoint);
    parent21=parent2(1:crossoverpoint);
    
    s1=[];
    s2=[];
    
    for j=1:numel(parent1)
       
        if sum(parent1(j)~=parent21)==numel(parent21);
            s2=[s2 parent1(j)];
        end
        
        if sum(parent2(j)~=parent12)==numel(parent12);
            s1=[s1 parent2(j)];
        end
        
    end
    
    offspring1=[parent12 s1];
    offspring2=[parent21 s2];
    
    
    % evaluate fitness
    
    offspring1(Ncities+1)=0;
    offspring2(Ncities+1)=0;
    
    for z=1:(Ncities-1)
        offspring1(Ncities+1)=offspring1(Ncities+1)+sqrt(((cities(offspring1(z+1),2)-cities(offspring1(z),2)).^2)+(cities(offspring1(z+1),3)-cities(offspring1(z),3)).^2);
    end
    
    for z=1:(Ncities-1)
       offspring2(Ncities+1)=offspring2(Ncities+1)+sqrt(((cities(offspring2(z+1),2)-cities(offspring2(z),2)).^2)+(cities(offspring2(z+1),3)-cities(offspring2(z),3)).^2);
    end        
             
    
    offspringset=[offspringset;offspring1;offspring2];
       
    end
    
    
    %%% end of crossover
    
    
    %% Mutation
    
     % selecting a parent for mutation
    m=randperm(popsize);
    
    for j=1:numberofmutation
        individualsForMutation(j,:)=pop(m(j),1:Ncities);
    end
    
   
    %doing the mutation
    Mutatedset=[];
    for z=1:numberofmutation
       
        parent=individualsForMutation(z,:);
        
        mutationPoint1=ceil(Ncities*rand);
        mutationPoint2=ceil(Ncities*rand);
        
        temp=parent(mutationPoint1);
        parent(mutationPoint1)=parent(mutationPoint2);
        parent(mutationPoint2)=temp;
        
        
        % Evaluate fitness of mutateed members
        
        parent(Ncities+1)=0;
        
        for z=1:(Ncities-1)
            parent(Ncities+1)=parent(Ncities+1)+sqrt(((cities(parent(z+1),2)-cities(parent(z),2)).^2)+(cities(parent(z+1),3)-cities(parent(z),3)).^2);
        end          
        
        Mutatedset=[Mutatedset;parent];
               
        
    end
%%% End mutation

%new population :

pop2=[selected;offspringset;Mutatedset];

        
   %sort the new population     
    for z=1:size(pop)
        for j=1:z
            if pop(j,Ncities+1)>pop(z,Ncities+1)
                temp=pop(z,:);
                pop2(z,:)=pop(j,:);
                pop2(j,:)=temp;
            end
        end
    end
    %find the best root
    
    for j=1:popsize
        pop(j,:)=pop2(j,:);
    end
    
    bestfit(i)=pop(1,Ncities+1);
    
    %calculate  mean
    mean=0;
    for j=1:popsize
       mean=mean+pop(j,Ncities+1);
    end
    meanfit(i)=mean/popsize;    
        
        
end   

best_root=pop(1,1:8)
min_distance=pop(1,9)
subplot(2,1,1)
plot(bestfit,'.r')

subplot(2,1,2)
plot(meanfit)
    

