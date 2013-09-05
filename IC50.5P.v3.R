
require(ggplot2)
require(stats)
op <- par(no.readonly = TRUE)

IC50.5P <- function(dose, Resp, T0 = NA, Ctrl = NA, LPweight = 0.25, fixB = NA, fixT = NA, fixS = NA,
#                    Plot = TRUE, pcol = 'grey50', lcol = 'grey25', unit="nM", Title="",
                    AddIC = .5, AddSd = TRUE, output = TRUE,...)
  {
  
  object <- cellResp(dose = dose, Resp = Resp, LPweight = LPweight)
  object@survProp <- .survProp(Resp, T0, Ctrl)
	# Optimisation step using nlm()
      init <- .initPar(.getDose(object), getSurvProp(object))
      model4 <- nlm(f = .sce.5P, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = 1)
      model4$estimate
      model5 <- nlm(f = .sce.5P, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = fixS)
      model5$estimate
  
	# Get best model
      bestModel <- .getBestModel(object, model4, model5)
      Param <- bestModel$param
      object@goodness <- bestModel$goodness
      object@model <- bestModel$model
  
	# Estimate critical points
		  object@xCurve <- seq(min(dose, na.rm = TRUE)*1.1, max(dose, na.rm = TRUE)*1.1, length = 100)
      object@yCurve <- .LP5(Param$bottom, Param$top, Param$xmid, Param$scal, Param$s, getXcurve(object))

	    # Inflexion point using the second derivative
		    Xflex = Param$xmid + (1/Param$scal)*log10(Param$s)
				Yflex = Param$bottom + (Param$top - Param$bottom)*(Param$s/(Param$s+1))^Param$s

			# Compute simulations to estimate the IC50 conf. interval
        P <- seq(.1, .9, by = .1)
        sigma <- summary(getGoodness(object))$sigma
        estimates <- lapply(P, function(p){.estimateRange(p, sigma, Param = Param, B = 1e4)})
        estimates <- cbind.data.frame(Resp = P, do.call(rbind, estimates))
        colnames(estimates) <- c('Surv', 'Dmin', 'D', 'Dmax')
        object@estimates <- estimates

# 	# Graphics
#     if(Plot){
#       par(las = 1, cex.axis = 1.5, cex.lab = 1.75, mar = c(6.5, 5.5, 4, 2), mgp = c(3.5, 1, 0))
#       plot(object, pch = 19)
#       par(op)
#     }
  return(object)
}



## Demo
source('~/Documents/MyProjects/FredScripts/LP5_classDef.R')
source('~/Documents/MyProjects/FredScripts/LP5_functions.R')
op <- par(no.readonly = TRUE)

x <- seq(log10(0.001), log10(10), len = 8)
y <- lapply(1:3, function(i) {rev(.LP5(0.1, 2.8, median(x), 1.5, 1, x) + rnorm(length(x), sd = .25))})
Resp <- do.call(c, y) + rnorm(3*8, .5, .1)
dose <- rep(signif(10^x, 2), 3)
dose <- log10(dose)

test <- IC50.5P(dose, Resp, T0 = .5, Ctrl = max(Resp))
plot(test)
getEstimates(test)


