genalg <-function(formula, data, family=gaussian(), population=200, generations=100, mutateRate=0.02, elitism=NA, zeroOneRatio=10, verbose=TRUE, fun.eval=evalAICc) {
    # Get the exponential family of the response
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
    # Match the variables that are referenced in the function call
    mf = match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    
    # We need to evaluate the function call in the next environment up because we've added a level of abstraction with this function.
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf)
    mt <- attr(mf, "terms")
    
    # Get the data and the weights
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    w <- model.weights(mf)
    
    #Figure out whether there's an intercept, and which columns of x are dummy variables for factors:
    intercept = ifelse(attr(mt, "intercept")==0, FALSE, TRUE)
    var.id = attr(x, "assign")
    
    # Check for problems with the (prior) weights
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    if (is.null(w)) 
        w <- rep(1, nrow(data))
    if (any(is.na(w))) 
        stop("NAs in weights")
    if (any(w < 0)) 
        stop("negative weights")
    
    #Create the object that will hold the output
    result = list()
    
    #Drop any rows with NA values
    na.rows = (which(is.na(data))-1) %% dim(data)[1] + 1
    if (length(na.rows)>0)
        data = data[-na.rows,]
        
    p = ncol(data)
    n = nrow(data)
    
    #Maximum number of predictor variables:
    result[['ga']] = rbga.bin(x, y, w, intercept, family, zeroToOneRatio=zeroOneRatio, evalFunc=fun.eval, monitorFunc=monitor, mutationChance=mutateRate, n.pop=population, n.iter=generations, verbose=verbose)
    
    indx = which.min(result[['ga']]$evaluations)
    indiv = as.logical(drop(result[["ga"]]$population[indx,]))
    
    #Assemble the design matrix for the highest-fitness model:
    cols = which(var.id %in% which(indiv==1))
    x.final = x[,cols]
    if (intercept) {x.final = cbind(1, x.final)}
    
    result[['vars']] = attr(mf, "term.labels")[which(indiv==1)]
    #result[['formula']] = as.formula(paste(response.name, "~", paste(result[['vars']], collapse="+"), sep=""))
    result[["model"]] = glm.fit(y=y, x=x.final, weights=w, family=family)
    
    result[['fitted']] = fitted(result[['model']])
    result[['residuals']] = residuals(result[['model']])
    result[['actual']] = result[['fitted']] + result[['residuals']]
    
    class(result) = "genalg"
    
    result
}
