function roulette_wheel(p)
    threshold=rand()*sum(p)
    distrib=cumsum(p)
    return findfirst(threshold.<=distrib)
end

function rank_selection(n)
    # The population is supposed sorted: last=highest fitness
    # n size of population
    p=convert(Array,1:n)
    i,j=roulette_wheel.([p, p])
    return i,j
end

function fitness_proportionate_selection(fitnesses)
    i,j=roulette_wheel.([fitnesses,fitnesses]) 
    return i,j
end

function selection_operator(fitnesses,proba)
    # fitnesses are sorted: lowest to highest
    operator=roulette_wheel(proba)
    if operator==1
        i,j=rank_selection(length(fitnesses)) 
    else
        i,j=fitness_proportionate_selection(fitnesses)
    end
    return i,j
end