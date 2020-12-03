function [mating_pool]=g_roulette_wheel(pop,n,P_select)
    cdf(1)=P_select(1);
	for itr=2:n
		cdf(itr)=cdf(itr-1)+P_select(itr);
	end

	for itr=1:n
		for j=1:n
			if rand<=cdf(j)
				mating_pool(itr,:)=pop(j,:);
				break;
            end
        end
    end
return;