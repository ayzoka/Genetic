function [best]=genetic(N,Pc,Pm,ITER,BS,L,cnt)

    
    m=3;
    Lo=[-4 -1.5 -1];
    Hi=[2 1 1];
    
    Population=round(rand(N,L));
    best_so_far=[];
    average_fitness=[];
    
    for i=1:ITER
          Pc=Pc-0.001;
          Pm=Pm+0.0005;
          Pc=max(Pc,0.49);
          Pm=min(Pm,0.1);
        [real_val]=chrom_decode(Population,N,L,BS,m,Lo,Hi);
        [selection_probability,~,ave_Fit,max_fit,opt_sol]=fit_eval(real_val,N,m);

        if(i==1)
            best_so_far(i)=max_fit;
            final_sol=opt_sol;
        elseif(max_fit>best_so_far(i-1))
            best_so_far(i)=max_fit;
            final_sol=opt_sol;
        else
            best_so_far(i)=best_so_far(i-1);
        end
        
        average_fitness(i)=ave_Fit;
        [mating_pool]=g_roulette_wheel(Population,N,selection_probability);
        [new_pop]=g_crossover(mating_pool,Pc,N,L);
        [Population]=g_mutation(new_pop,Pm,N,L);
    end

    
    best=best_so_far(ITER);
    disp("result: "+num2str(best));
    disp(final_sol);
    x=1:ITER;
    f=figure('visible','off');
    plot(x,best_so_far,'k',x,average_fitness(1:2000),'.-k');
    xlabel('generation');
    ylabel('fitness function')
    legend('best so far','average fitness')
    saveas(f,"multimodalbinarycoded"+cnt,'svg')

return;

