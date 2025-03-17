library(ReIns)
library(actuar)

MethodNormal <- function(lambda, alpha, xm) {
  moments <- GetParetoMoments(alpha, xm)
  aggmoments <- GetMomentsOfAggrDistr(lambda, lambda, moments[1], moments[2])
  Fs <- aggregateDist("normal", moments = aggmoments)
return(c(mean(Fs), quantile(Fs, 0.95)))
}

MethodConvolve <- function(lambda, alpha, xm) {
  print("convolve")
  xlim <- lambda*10
  freq <- discretize(ppois(x, lambda = lambda), from = 0, to = xlim, by = max(round(xlim/20), 1))
  #freq <- freq[freq > 0.0001]
  freq <- c(1-sum(freq), freq)
  sev <- discretize(ReIns::ppareto(x, shape = alpha, scale = xm), from = 0, to = 10*(xm*5), by = xm/100)
  sev <- c(1-sum(sev), sev)
  Fs <- aggregateDist("convolution", model.freq = freq, model.sev = sev, x.scale = xm/100)
return(c(mean(Fs), quantile(Fs, 0.95)))
}

MethodRecursive <- function(lambda, alpha, xm) {
  print("recurs")
  freqred <- (max(ceiling(log(lambda/16)/log(2)), 0))
  print(freqred)
  sev <- discretize(ReIns::ppareto(x, shape = alpha, scale = xm), from = 0, to = 10*(xm*20), by = xm/100)
  sev <- c(1-sum(sev), sev)
  Fs <- aggregateDist("recursive", model.freq = "poisson", model.sev = sev, lambda = lambda/2^{freqred}, x.scale = xm/100, maxit = 100000, convolve = freqred)  
return(c(mean(Fs), quantile(Fs, 0.95)))
}

MethodSimulation <- function(lambda, alpha, xm) {
  print("sim")
  lambda <<- lambda
  alpha <<- alpha
  xm <<- xm
  model.freq <- expression(data = rpois(lambda))
  model.sev <- expression(data = ReIns::rpareto(alpha, xm))
  numiter <- round(GetNumberOfIterationsFull(lambda, alpha, xm))
    
  Fs <- aggregateDist("simulation", nb.simul = numiter, model.freq, model.sev)
return(c(mean(Fs), quantile(Fs, 0.95), numiter))
}

MethodSimulationGlob <- function(lambda, alpha, xm) {
  print("method global")
  lambda <<- lambda
  alpha <<- alpha
  xm <<- xm
  model.freq <- expression(data = rpois(lambda))
  model.sev <- expression(data = ReIns::rpareto(alpha, xm))
  numiter <- round(20000000/max(lambda, 1))
  print(numiter)  
  Fs <- aggregateDist("simulation", nb.simul = numiter, model.freq, model.sev)
return(c(mean(Fs), quantile(Fs, 0.95), numiter))
}


GetParetoMoments <- function(alpha, xm) {
  if (alpha <= 1) {
    mu <- Inf
  } else {
    mu <- alpha*xm/(alpha-1)
  }
  if (alpha <= 2) {
    sigma <- Inf
  } else {
    sigma <- xm^2*alpha/((alpha-2)*(alpha-1)^2)
  #  sigma <- (alpha*xm^2/(alpha-2))
  }
c(mu, sigma)
}

GetMomentsOfAggrDistr <- function(meanfreq, varfreq, meansev, varsev) {
c(meanfreq*meansev, meanfreq*varsev + varfreq*(meansev)^2)
}

GetMoments <- function(lambda, alpha, xm) {
  moments <- GetParetoMoments(alpha, xm)
  GetMomentsOfAggrDistr(lambda, lambda, moments[1], moments[2])
}

GetNumberOfIterations <- function(mu, sigma) {
  sigma*((2/(0.01*mu))^2)
}

GetNumberOfIterationsFull <- function(lambda, alpha, x0, lowbignum = 5000, topbignum = 20000000) {
  moments <- GetParetoMoments(alpha, x0)
  aggmoments <- GetMomentsOfAggrDistr(lambda, lambda, moments[1], moments[2])
  BigNum <- GetNumberOfIterations(aggmoments[1], aggmoments[2])
  max(min(BigNum, round(topbignum/max(lambda, 1))), lowbignum)
}

StudyOneCase <- function(lambda, alpha, xm) {
  timeconv <- system.time(convres <- MethodConvolve(lambda, alpha, xm))
  gc()
  timerec <- system.time(recres <- MethodRecursive(lambda, alpha, xm))
  gc()
  timesim <- system.time(simres <- MethodSimulation(lambda, alpha, xm))
  gc()
  timeglob <- system.time(globres <- MethodSimulationGlob(lambda, alpha, xm))
  gc()
  moments <- GetParetoMoments(alpha, xm)
  aggmoments <- GetMomentsOfAggrDistr(lambda, lambda, moments[1], moments[2])
  data.frame(timeconv = timeconv[3],
             timerec = timerec[3],
			 timesim = timesim[3],
			 timeglob = timeglob[3],
			 truemean = aggmoments[1],
			 convresmean = convres[1],
			 convresvar = convres[2],
			 recresmean = recres[1],
			 recresvar = recres[2],
			 simresmean = simres[1],
			 simresvar = simres[2],
			 globresmean = globres[1],
			 globresvar = globres[2])
}

StudyAllCases <- function() {
  lambdas <- c(0.1, 1, 2, 5, 10, 15, 30, 50, 75, 100)
  alphas <- c(2.5, 4, 10)
  xms <- c(100000, 1000000, 10000000)
  len <- length(alphas)*length(xms)*length(lambdas)
  resframe <- data.frame(timeconv = numeric(len),
                         timerec = numeric(len),
						 timesim = numeric(len),
						 timeglob = numeric(len),
						 truemean = numeric(len),
						 convresmean = numeric(len),
						 convresvar = numeric(len),
						 recresmean = numeric(len),
						 recresvar = numeric(len),
						 simresmean = numeric(len),
						 simresvar = numeric(len),
						 globresmean = numeric(len),
						 globresvar = numeric(len))
  id <- 1					 
  for (lambda in lambdas) {
    for (alpha in alphas) {
	  for (xm in xms) {
	    print(c(lambda, alpha, xm))
	    resframe[id, ] <- StudyOneCase(lambda, alpha, xm)
		id <- id + 1
		write.table(resframe, "C:\\RScripts\\res.csv", sep = ';', row.names = FALSE)
		gc()
	  }
	}
  }
resframe  
}

StudyNormal <- function() {
  lambdas <- c(0.1, 1, 2, 5, 10, 15, 30, 50, 75, 100)
  alphas <- c(2.5, 4, 10)
  xms <- c(100000, 1000000, 10000000)
  len <- length(alphas)*length(xms)*length(lambdas)
  resframe <- data.frame(timenorm = numeric(len),
						 normmean = numeric(len),
						 normvar = numeric(len))
  id <- 1					 
  for (lambda in lambdas) {
    for (alpha in alphas) {
	  for (xm in xms) {
	    print(c(lambda, alpha, xm))
	    resframe[id, ] <- StudyOneCaseNormal(lambda, alpha, xm)
		id <- id + 1
		write.table(resframe, "C:\\RScripts\\res2.csv", sep = ';', row.names = FALSE)
		gc()
	  }
	}
  }
resframe  
}

StudyOneCaseNormal <- function(lambda, alpha, xm) {
  timenormal <- system.time(normres <- MethodNormal(lambda, alpha, xm))
  gc()
  data.frame(timenormal = timenormal[3],
			 normmean = normres[1],
			 normvar = normres[2])
}
