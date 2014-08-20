evalAICc <-
function(model, tol=1e-4) {
    returnVal = Inf

    #if (model$deviance > tol) {
        n = length(model$residuals)
        df = model$rank
        returnVal = model$aic + 2*df*(df+1)/(n-df-1)
    #}
    
    if (is.na(returnVal)) { returnVal=Inf }
    
    return(returnVal)
}
