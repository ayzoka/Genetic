function [mating_pool]=g_roulette_wheel(Population,N,selection_probablity)
	cdf(1)=selection_probablity(1);
	for i=2:N
		cdf(i)=cdf(i-1)+selection_probablity(i);
	end

	for i=1:N
		q=rand();
		for j=1:N
			if q<=cdf(j)
				mating_pool(i,:)=Population(j,:);
				break;
			end 
		end
	end
	return;