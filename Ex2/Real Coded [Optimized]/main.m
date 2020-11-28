clc;
clear;
close all;
m=3;
lwb=[0 14.7 0];
upb=[15.1 94.2 5371];

iter=500;
population_size=50;

population=zeros(population_size,m);
temp_population=zeros(population_size,m);
mu=0.0015;

sorted_population= zeros (population_size,m+1);
F=zeros(population_size,1);
cnt_crossover=4;
cnt_mutation=5;

 for i=1:population_size
    for n=1:m
        population(i,n)= unifrnd(lwb(1,n),upb(1,n));
    end
    F(i)= F_eval(population(i,:));
 end
 
sorted_population(:,1:m)= population(:,:);
sorted_population(:,m+1)=F;
sorted_population= sortrows(sorted_population,m+1);
 
% Generation
for ii=1:1:iter
    k=1;

    temp_population(k,1:m)=sorted_population(1,1:m);
    k=k+1;

    for j=1:cnt_crossover
        y1(j) = geornd(0.1)+1;
        while   y1(j)> population_size
                y1(j) = geornd(0.1)+1;
        end
        y2(j) = geornd(0.1)+1;
        while   y2(j)> population_size
                y2(j) = geornd(0.1)+1;
        end
    end
    for u=1:cnt_crossover
        parent1= sorted_population(y1(u),1:m);
        parent2= sorted_population(y2(u),1:m);

        [Children]= arithmetic_crossover (parent1, parent2);
        temp_population(k,1:m)=Children(1,:);

        temp_population(k,1:m)=max(temp_population(k,1:m),lwb);
        temp_population(k,1:m)=min(temp_population(k,1:m),upb);

        k=k+1;
        temp_population(k,1:m)=Children(2,:);

        temp_population(k,1:m)=max(temp_population(k,1:m),lwb);
        temp_population(k,1:m)=min(temp_population(k,1:m),upb);

        k=k+1;

    end
    
    for e=1:cnt_mutation
        parent=sorted_population(unidrnd(population_size),1:m);
        [child]= g_mutation (parent,mu,lwb,upb);

        temp_population(k,1:m)=child;
        temp_population(k,1:m)=max(temp_population(k,1:m),lwb);
        temp_population(k,1:m)=min(temp_population(k,1:m),upb);

        k=k+1;
    end

    for k=k:1:population_size
        replicated_child=sorted_population(randi([1 population_size]),1:m);
        temp_population(k,1:m)=replicated_child;
        k=k+1;
    end

    for iii=1:1:population_size
        F(iii)= F_eval(temp_population(iii,:));
    end

    sorted_population(:,1:m)= temp_population;
    sorted_population(:,m+1)=F(:,:);
    sorted_population= sortrows(sorted_population,m+1);

    best_res(ii)= sorted_population(1,m+1);
end
plot(best_res,'b');
pause(0.1)
xlabel('Generation')
ylabel('Best Value')
grid on
Xmin=sorted_population(1,1:m)
Fval=best_res(1,iter)
