#' Run the genetic algorithm.
#' 
#' @param x model.matrix 
#' @param y response vector
#' @param w weight vector
#' @param intercept use an intercept?
#' @param family \code{family} object to which the response belongs
#' @param n.pop number of individuals in each generation
#' @param n.iter number of generations
#' @param mutationChance probability of a mutation at each location on a chromosome (when that chromosome is eligible for mutation)
#' @param elitism number of highest-fitness individuals from each generation to pass on directly to the next generation, with no chance for mutation
#' @param zeroToOneRatio odds of getting a 1 at each location on a chromosome that is being randomly generated, either from mutation or from initial generation
#' @param monitorFunc function that monitors progress in the algorithm and prints a summary after each generation
#' @param fitness function that evaulates the fitness for a given chromosome 
#' @param showSettings print the settings at the start of the genetic algorithm? default is \code{FALSE}.
#' @param verbose print detailed information about the progress of the genetic algorithm? Default is \code{FALSE}.
#' 
rbga.bin = function(x, y, w, intercept, family, suggestions=NULL, n.pop=200, n.iter=100, mutationChance=NA, elitism=NA, zeroToOneRatio=10, monitorFunc=NULL, fitness=NULL, showSettings=FALSE, verbose=FALSE) {
    if (is.na(mutationChance)) {mutationChance = 1/(n.var + 1)}
    if (is.na(elitism)) {elitism = floor(n.pop/5)}
    
    #Check some conditions that would abort this run:
    if (verbose) {cat("Testing the sanity of parameters...\n")}
    if (is.null(fitness)) {stop("A evaluation function must be provided. See the fitness parameter.")}
    if (n.pop < 5) {stop("The population size must be at least 5.")}
    if (n.iter < 1) {stop("The number of iterations must be at least 1.")}
    if (!(elitism < n.pop)) {stop("The population size must be greater than the elitism.")}
    
    #Intercept is an always-or-never proposition:
    var.id = attr(x, "assign")
    n.var = length(unique(var.id))
    if (intercept) {n.var = n.var - 1}
    
    if (showSettings) {
        if (verbose) {cat("The start conditions:\n")}
        result = list(size=n.var, suggestions=suggestions, 
                      n.pop=n.pop, n.iter=n.iter, elitism=elitism, 
                      mutationChance=mutationChance)
        class(result) = "rbga"
        cat(summary(result))
    } else {
        if (verbose) {cat("Not showing GA settings...\n")}
    }
    
    if (n.var > 0) {
        if (!is.null(suggestions)) {
            if (verbose) {cat("Adding suggestions to first population...\n")}
            population = matrix(nrow=n.pop, ncol=n.var)
            suggestionCount = dim(suggestions)[1]
            for (i in 1:suggestionCount) { population[i,] = suggestions[i,] }
            
            if (verbose) cat("Filling others with random values in the given domains...\n")
            for (child in (suggestionCount + 1):n.pop) {
                population[child, ] = rbinom(n.var, p=1/(1+zeroToOneRatio), size=1)
                while (sum(population[child, ]) == 0) {
                    population[child, ] = rbinom(n.var, p=1/(1+zeroToOneRatio), size=1)
                }
            }
        } else {
            if (verbose) {cat("Starting with random values in the given domains...\n")}
            population = matrix(nrow=n.pop, ncol=n.var)
            
            #Produce a complete generation:
            for (child in 1:n.pop) {
                population[child,] = rbinom(n.var, p=1/(1+zeroToOneRatio), size=1)
                while (sum(population[child,]) == 0) {
                    population[child,] = rbinom(n.var, p=1/(1+zeroToOneRatio), size=1)
                }
            }
        }
        
        #Set up variables to hold results:
        bestEvals = rep(NA, n.iter)
        meanEvals = rep(NA, n.iter)
        evalVals = rep(NA, n.pop)
        
        #Loop over the generations:
        for (iter in 1:n.iter) {
            if (verbose) {cat(paste("Starting iteration", iter, "\n"))}
            if (verbose) {cat("Calucating evaluation values... ")}
            
            #Evaluate all individuals in this generation:
            for (j in 1:n.pop) {
                #Evaluate an individual model:
                if (is.na(evalVals[j])) {
                    cols = which(var.id %in% which(population[j,]==1))    
                    x.indiv = x[,cols]
                    if (intercept) {x.indiv = cbind(1, x.indiv)}
                    
                    model = glm.fit(x=x.indiv, y=y, weights=w, family=family)
                    evalVals[j] = fitness(model)
                    if (verbose) {cat(".")}
                }
            }
            
            #Summarize the evaluated models:
            bestEvals[iter] = min(evalVals)
            meanEvals[iter] = mean(evalVals)
            if (verbose) {cat(" done.\n")}
            
            
            if (!is.null(monitorFunc)) {
                if (verbose) {cat("Sending current state to rgba.monitor()...\n")}
                result = list(type="binary chromosome", size=n.var, 
                              n.pop=n.pop, iter=iter, n.iter=n.iter, 
                              population=population, elitism=elitism, 
                              mutationChance=mutationChance, evaluations=evalVals, 
                              best=bestEvals, mean=meanEvals)
                class(result) = "rbga"
                monitorFunc(result, output=output, data=data)
            }
            
            #If this isn't the last generation, then prime the next one from the current generation
            if (iter < n.iter) {
                if (verbose) {cat("Creating next generation...\n")}
                newPopulation = matrix(nrow = n.pop, ncol = n.var)
                newEvalVals = rep(NA, n.pop)
                
                #Sort individuals based on their fitness
                if (verbose) {cat("  sorting results...\n")}
                sortedEvaluations = sort(evalVals, index = TRUE)
                sortedPopulation = matrix(population[sortedEvaluations$ix,], ncol=n.var)
                
                #Elitism means copying the most fit individuals directly into the new generation:
                if (elitism > 0) {
                    if (verbose) {cat("  applying elitism...\n")}
                    newPopulation[1:elitism,] = sortedPopulation[1:elitism,]
                    newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
                }
                
                #Once past the "elite", we have random(?) mating:
                if (n.var > 1) {
                    if (verbose) {cat("  applying crossover...\n")}
                    
                    
                    for (child in (elitism + 1):n.pop) {
                        #Draw two parents for this child
                        #parentProb = dnorm(1:n.pop, mean=0, sd=(n.pop/3))
                        #parentIDs = sample(1:n.pop, 2, prob=parentProb)
                        parentIDs = sample(1:n.pop, 2)
                        parents = sortedPopulation[parentIDs,]
                        
                        #Crossover: from 1:crossover, copy the first parent's chromosome.
                        #From (crossover+1):n.var, copy the second parent's chrmosome.
                        crossOverPoint = sample(0:n.var, 1)
                        
                        if (crossOverPoint == 0) {
                            #In this case, we're copying the second parent exactly.
                            #We already know its fitness, too.
                            newPopulation[child,] = parents[2,]
                            newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
                        } else if (crossOverPoint == n.var) {
                            #Copy the first parent exactly. We know its fitness already.
                            newPopulation[child,] = parents[1,]
                            newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
                        } else {
                            #Some combination of the parents:
                            newPopulation[child,] = c(parents[1,][1:crossOverPoint], parents[2,][(crossOverPoint+1):n.var])
                            
                            #If we produced a totally blank chromosome, then generate a new one randomly
                            while (sum(newPopulation[child, ]) == 0) {
                                newPopulation[child,] = rbinom(n.var, p=1/(1+zeroToOneRatio), size=1)
                            }
                        }
                    }
                } else {
                    if (verbose) {cat("  cannot crossover (#n.var=1), using new randoms...\n")}
                    newPopulation[(elitism+1):n.pop,] = sortedPopulation[sample(1:n.pop, n.pop-elitism),]
                }
                
                population = newPopulation
                evalVals = newEvalVals
                
                #Now maybe mutate (the elites avoid mutation):
                if (mutationChance > 0) {
                    if (verbose) {cat("  applying mutations... ")}
                    mutationCount = 0
                    
                    #Find mutations for each individual:
                    for (individual in (elitism + 1):n.pop) {
                        mutate = rbinom(n.var, p=mutationChance, size=1)
                        n.mutate = sum(mutate)
                        population[individual, which(mutate==1)] = rbinom(n.mutate, p=1/(1+zeroToOneRatio), size=1)
                        mutationCount = mutationCount + n.mutate
                    }
                    if (verbose) {cat(paste(mutationCount, "mutations applied\n"))}
                }
            }
        }
    }
    
    result = list(type="binary chromosome", size=n.var, n.pop=n.pop, 
                  n.iter=n.iter, suggestions=suggestions, population=population, 
                  elitism=elitism, mutationChance=mutationChance, evaluations=evalVals, 
                  best=bestEvals, mean=meanEvals)
    class(result) = "rbga"
    return(result)
}
