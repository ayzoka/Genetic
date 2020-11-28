function [Children]= arithmetic_crossover (p1, p2)
    alpha=rand;
    c1=alpha*p1+(1-alpha)*p2;
    c2=(1-alpha)*p1+alpha*p2;
    Children=[c1;c2];
end
