genetic<-function(formula, data, family=gaussian, weights, selection.criterion, population=200, generations=100, mutateRate=0.02, zeroOneRatio=10, verbose=TRUE, na.action=NULL, subset=NULL, offset=NULL, ...) {
    #Interpret the provided arguments:
    call = match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    
    #Evaluate the family:
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
    #Create the object that will hold the output
    result = list()
    result[['call']] = call
    result[['formula']] = as.formula(formula, env=data)
    result[['selectonly']] = selectonly
    result[['selection.criterion']] = selection.criterion = match.arg(selection.criterion)
    
    #Get the model data:
    Y = model.response(mf)
    X = model.matrix(mf, data=mf, contrasts)[,-attr(mt, "intercept")]
    vars = attr(mt, "term.labels")
    
    #Get the observation weights
    if (missing(weights)) w = rep(1, nrow(X))
    else w = as.vector(model.weights(mf))
    
    #Create the object that will hold the output
    result = list()
    
    #Drop any rows with NA values
    na.rows = (which(is.na(data))-1) %% dim(data)[1] + 1
    if (length(na.rows)>0)
        data = data[-na.rows,]
        
    m = ncol(X)
    n = nrow(X)
    
    #Maximum number of predictor variables:
    m = ncol(data) - 1
    result[['ga']] = rbga.bin(X, Y, size=m, zeroToOneRatio=zeroOneRatio, evalFunc=evalAICc, monitorFunc=monitor, mutationChance=mutateRate, popSize=population, iters=generations, verbose=verbose)
    
    indx = which.min(result[['ga']]$evaluations)
    indiv = as.logical(drop(result[["ga"]]$population[indx,]))
    
    result[['vars']] = predictor.names[indiv]
    result[['formula']] = as.formula(paste(response.name, "~", paste(result[['vars']], collapse="+"), sep=""))
    result[["model"]] = lm(formula=result[['formula']], data=data)
    
    result[['fitted']] = fitted(result[['model']])
    result[['residuals']] = residuals(result[['model']])
    result[['actual']] = result[['fitted']] + result[['residuals']]
    
    class(result) = "genetic"
    result
}
