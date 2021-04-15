include("Selection Operators.jl")

function horizontal_split(parent1,parent2)
    split_at=rand(2:size(parent1,1))
    offspring1=copy(parent1)
    offspring2=copy(parent2)
    offspring1[1:split_at-1,:]=parent2[1:split_at-1,:]
    offspring2[1:split_at-1,:]=parent1[1:split_at-1,:]
    return offspring1,offspring2
end

function vertical_split(parent1,parent2)
    split_at=rand(2:size(parent1,1))
    offspring1=copy(parent1)
    offspring2=copy(parent2)
    offspring1[:,1:split_at-1]=parent2[:,1:split_at-1]
    offspring2[:,1:split_at-1]=parent1[:,1:split_at-1]
    return offspring1,offspring2    
end

function uniform_crossover(parent1,parent2)
    λ=rand()
    offspring1=λ*parent1+(1-λ)*parent2
    offspring2=λ*parent2+(1-λ)*parent1
    return offspring1,offspring2    
end

crossover_operators=[uniform_crossover, vertical_split, horizontal_split]

function crossover_operator(parent1,parent2,rate,proba)
    if rand()<rate
        i=roulette_wheel(proba)
        parent1,parent2=crossover_operators[i](parent1,parent2)
    end
    return parent1,parent2
end