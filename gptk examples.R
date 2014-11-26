library(gptk)
demOptimiseGp

####this is the guts of demOptimiseGp
path = getwd()
filename = "demOptimiseGp"
png = FALSE 
gif = FALSE
  options = gpOptions()
  options$kern$comp = list("rbf", "white")
  lengthScale = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
  figNo = 0
  x = matrix(seq(-1, 1, length = 6), ncol = 1)
  xtest = matrix(seq(-1.5, 1.5, length = 200), ncol = 1)
  
  trueKern = kernCreate(x, list(type = "cmpnd", comp = list("rbf", 
                                                            "white")))
  kern = trueKern
  K = kernCompute(trueKern, x)
  y = t(gaussSamp(matrix(0, 6, 1), K, 1))
  y = scale(y, scale = FALSE)
  model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)
  graphics.off()
  ll = c()
  llLogDet = c()
  llFit = c()
  for (i in 1:length(lengthScale)) {
    inithypers = log(c(1/(lengthScale[i]), 1, model$kern$comp[[2]]$variance))
    model = gpExpandParam(model, inithypers)
    invK = model$invK_uu
    logDetK = model$logDetK_uu
    ll[i] = gpLogLikelihood(model)
    llLogDet[i] = -0.5 * (logDetK + dim(y)[1] * log(2 * pi))
    llFit[i] = -0.5 * t(y) %*% invK %*% y
    meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return = TRUE)
    dev.new()
    plot.new()
    gpPlot(model, xtest, meanVar$mu, meanVar$varsigma, ylim = c(-2.5, 
                                                                2.5), xlim = range(xtest), col = "black")
    figNo = figNo + 1
    if (png) {
      dev.copy2eps(file = paste(path, "/", filename, "1_", 
                                as.character(figNo), ".eps", sep = ""))
      system(paste("eps2png ", path, "/", filename, "1_", 
                   as.character(figNo), ".eps", sep = ""))
    }
    dev.new()
    plot.new()
    matplot(lengthScale[1:i], cbind(ll[1:i], llLogDet[1:i], 
                                    llFit[1:i]), type = "l", lty = c(1, 2, 3), log = "x", 
            xlab = "length-scale", ylab = "log-probability")
    legend(x = "topleft", c("marginal likelihood", "minus complexity penalty", 
                            "data-fit"), lty = c(1, 2, 3), col = c("black", "red", 
                                                                   "green"))
    if (png) {
      dev.copy2eps(file = paste(path, "/", filename, "2_", 
                                as.character(figNo), ".eps", sep = ""))
      system(paste("eps2png ", path, "/", filename, "2_", 
                   as.character(figNo), ".eps", sep = ""))
    }
  }
  if (gif) {
    system(paste("convert -delay 80 ", path, "/", filename, 
                 "1_*.png ", path, "/", filename, "1.gif", sep = ""))
    system(paste("convert -delay 80 ", path, "/", filename, 
                 "2_*.png ", path, "/", filename, "2.gif", sep = ""))
  }

