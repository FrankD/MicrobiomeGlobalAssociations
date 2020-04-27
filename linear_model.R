# Fitting simple linear marginal models to each species to determine the effect 
# size.
# We assume small p, and so the linear models can be fitted with lm().

# Main fitting function, one linear model per species.
fit_linear_models <- function(z, x) {
  S = dim(z)[2]
  
  a_sgn = matrix(NA, S, dim(x)[2])
  
  for(s in 1:S) {
    df = cbind(data.frame(z=z[,s]),
               x)

    linear_model = lm(z ~ ., data=df)
    
    a_sgn[s,] = sign(summary(linear_model)$coefficients[,1])
  }
  
  return(a_sgn)
}