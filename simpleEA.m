function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
%recorders
fitness_gen=[]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code

% Initialise a population
%% TODO
popsize = 4;
lowbound = 0;
highbound = 31;
pop = dec2bin(randi([lowbound,highbound],[popsize,1])-1);

% Evaluate the initial population
%% TODO
fitness = objective(bin2dec(pop));
[~,index] = max(fitness);
fitness_gen(1) = fitness(index);
solution_gen(1) = bin2dec(pop(index,:));
fitness_pop(1) = fitness(index);
nbGen = nbGen + 1;
nbEval = nbEval + popsize;

% Start the loop
while (nbEval<T)
% Reproduction (selection, crossver)
%% TODO
    % select 
    prob = fitness / sum(fitness);
    parent_index = randsrc(popsize/2,2,[1:popsize;prob']);
    % crossover
    new_pop = [];
    for k = 1:popsize/2 % only consider the popsize is even
        pc = randi(5);
        parent_1 = pop(parent_index(k,1),:);
        parent_2 = pop(parent_index(k,2),:);
        new_pop = [new_pop; [parent_1(1:pc-1),parent_2(pc:end)]];
        new_pop = [new_pop; [parent_2(1:pc-1),parent_1(pc:end)]];
    end
    
% Mutation
%% TODO
    for k=1:popsize
        if rand < 0.2
            p_index = randi(popsize);
            pm = randi(5);
            new_pop(p_index,pm) = dec2bin(~bin2dec(new_pop(p_index,pm)));
        end
    end
    
    nbGen = nbGen + 1;
    nbEval = nbEval + popsize;

    pop = new_pop;
    fitness = objective(bin2dec(pop));
    [~,index] = max(fitness);
    solution_gen = [solution_gen, bin2dec(pop(index,:))];
    fitness_pop = [fitness_pop, fitness(index)];
    fitness_gen = [fitness_gen, max(fitness(index), fitness_gen(nbGen-1))];

end

figure
plot(1:length(solution_gen),solution_gen)
title("Best Solution of Each Generation")

figure
plot(1:length(fitness_pop), fitness_pop)
title("Best Fitness in Each Generation")

figure
plot(1:length(fitness_gen), fitness_gen)
title("Best Fitness So Far")

[bestSoFarFit,best_index] = max(fitness_gen);
bestSoFarSolution = solution_gen(best_index);





