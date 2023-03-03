#-------------------------------------------------------------------------------
#                           Cahill QFC March 2023
# Questions:
#
# 1. Can we divine critter behavior as a latent state using  
# a distance correlated random walk?
# i.e., a fancy state-space spatial-temporal model
# 
# 2. Can we fit it to something with GLATOS-like data?
# 
# Some simulation play to start things off...
#-------------------------------------------------------------------------------

library(TMB)
library(mvtnorm) 
library(plotrix)

# function based on Auger-Methe et al. 2016
sim <- function(n, gamm, sdslo, sdsla, sdolo, sdola) {
  # n: length of track
  # gamm: autocorrelation/behavioral parameter, vector of length n
  # gamm control tortuosity
  # sdslo, sdsla: stand. dev. for process stochastiticy (randomness/unexplained variation in movement)
  # sdolo, sdola: stand. dev. in error term for measurements

  # Check that gamma is the length of the time series
  if (length(gamm) < n) gamm <- rep(gamm, length.out = n)
  x <- array(dim = c(n, 2))
  y <- array(dim = c(n, 2))

  # Variance-covariance matrix for the process error
  Sigma <- matrix(c(sdslo^2, 0, 0, sdsla^2), ncol = 2, byrow = TRUE)

  # Variance-covariance matrix for the measurment error
  SigmaO <- matrix(c(sdolo^2, 0, 0, sdola^2), ncol = 2, byrow = TRUE)

  # Simulate the first two true locations
  # i.e., condition on first two locations
  x[1, ] <- rmvnorm(1, c(0, 0), Sigma)
  x[2, ] <- rmvnorm(1, x[1, ], Sigma)

  # Simulate all other true locations
  for (i in 3:n) {
    # gamma changes the correlation of the random walk, controls both the speed and direction.
    x[i, ] <- rmvnorm(1, x[i - 1, ] + gamm[i] * (x[i - 1, ] - x[i - 2, ]), Sigma)
  }
  # Add measurement error to all true locs
  for (i in 1:n) {
    y[i, ] <- rmvnorm(1, x[i, ], SigmaO)
  }
  return(list(y = y, x = x))
}

# set some simulation parameter values:
sd_lat_p <- 0.1
sd_lon_p <- 0.1

nobs <- 20000 # total number of observations
behL <- nobs/2 # Number of steps to go through one behaviour

# Cycle of gamma values for each step.
gamm <- rep((cos(seq(-pi, pi, length.out = behL)) + 1.1) / 2.2, length.out = nobs)
plot(gamm, main = "Behaviorial state as per Auger-Methe et al. 2016")

# simulate some fake data
set.seed(100)
simdat <- sim(
  n = nobs, gamm = gamm, sdslo = sd_lon_p, sdsla = sd_lat_p,
  sdolo = .05, sdola = .05
)

# Blow holes in the simulated track by sampling it:

# how far apart are receivers (km)
grid_size <- 10

# calulate x and y range
x_range <- range(simdat$x[, 1])
y_range <- range(simdat$x[, 2])

# create a grid
x_seq <- seq(from = x_range[1], to = x_range[2], length.out = floor(diff(x_range) / grid_size))
y_seq <- seq(from = y_range[1], to = y_range[2], length.out = floor(diff(y_range) / grid_size))
receiver_locs <- expand.grid(x_seq, y_seq)

# calculate the detection probability for each true location | grid
beta1 <- 0.0025 # parameters represent moderate detection from Kraus et al.
beta2 <- 800
D <- 1:2000
det_fn <- 1 - (1 / (1 + 10^(-beta1 * (D - beta2))))
Dkm <- D / 1000
plot(det_fn ~ Dkm)

# function to generate detection probability given nearest receiver:
get_phat <- function(x){
  phat <- rep(NA, nrow(x))
  for(i in 1:nrow(simdat$x)){
    my_loc <- x[i,]
    dist <- min(dist(rbind(my_loc, receiver_locs)))
    phat[i] <- 1 - (1 / (1 + 10^(-beta1 * (dist*1000 - beta2))))
  }
  phat
}
phats <- get_phat(simdat$x)

# was a given location detected? 
idx <- rbinom(n = length(phats), size = 1, prob = phats)

# subset the simulated data based on what was detected:
simdat2 <- simdat$y[which(idx == 1),]
dt <- diff(which(idx == 1)) # calculate all the dt values

#-------------------------------------------------------------------------------
# estimate all of that in TMB...

# compile the TMB file
compile("dcrw/dcrw.cpp")

# create dynamically linked library:
dyn.load(dynlib("dcrw/dcrw"))

# set up data object
data <- list(y = t(simdat2), dt = dt)

# set up parameters to estimate--three fixed effects and a 
# random effect vector for gamma that is length(dt)
parameters <- list(logSdlat = 0, logSdlon = 0, logSdgamma = 0, 
                   gamma = rep(0, dim(data$y)[2]))

# create an automatic differentiation function, treat gamma as latent random effect
obj <- MakeADFun(data, parameters, random = "gamma", DLL = "dcrw")

# run the model
opt <- nlminb(obj$par, obj$fn, obj$gr)

# look at it
opt

# some checks
if(opt$convergence != 0 || any(abs(obj$gr(opt$par)) > 0.001)){message("model likely not converged!")}

# is hessian positive definite? 
sr <- summary(sdreport(obj)) 

# look at first few values
sr[1:10,] 

# extract some stuff
sts <- t(obj$env$parList()$gamma)

# Create a vector with estimated behavioural parameter gamma
gammEst <- sr[rownames(sr) %in% "gamma", 1]

# Look at the estimated parameters in comparison to simulated values
cbind(sr[c("sdLat", "sdLon"), ], SimValue = c(sd_lat_p, sd_lon_p))

sr[c("sdGamma"), ] # are standard errors too large?

#-------------------------------------------------------------------------------

# Plot the results
par(mfrow = c(1, 1))
plot(gammEst,
  pch = 19, cex = 0.15, ylim = c(-1, 1), xlab = "time",
  ylab = substitute(gamma[i]), 
  main = "gamma estimates + 95% CIs"
)
lines(gammEst + qnorm(0.975) * sr[rownames(sr) == "gamma", 2], 
      col = "steelblue", lwd = 0.3)
lines(gammEst + qnorm(0.025) * sr[rownames(sr) == "gamma", 2], 
      col = "steelblue", lwd = 0.3)

# Normalise the gamma estimate to get color intensity values between 0 and 1.
gammEstCol <- (gammEst - min(gammEst)) / (max(gammEst) - min(gammEst))

# Plot the two movement tracks side-by-side
png("dcrw/fig.png",
 width = 8, height = 4, units = "in", res = 2000
)
layout(matrix(1:2, nrow = 1))
# Simulated movement track
plot(simdat$x,
  cex = 0.35, pch = 19, col = rgb(gamm, 0, 0), axes = TRUE, xlab = "km", ylab = "km",
  main = expression(simulated ~ gamma[t] ~ and ~ receiver ~ grid), ty = "o"
) 

for (i in 1:nrow(receiver_locs)) {
  draw.circle(receiver_locs$Var1[i], receiver_locs$Var2[i],
              radius = 1.5,
              col = "white"
  )
}
points(receiver_locs, pch = 19, cex = 0.01, col = "black")

# Movement track with estimated gamma
plot(simdat$x,
  pch = 19, ty = "o", cex = 0.15, col = "grey", axes = TRUE, xlab = "km", ylab = "km",
  main = expression(estimated ~ gamma[t])
)
points(simdat2, pch = 19, cex = 0.35, col = rgb(gammEstCol, 0, 0)) 

dev.off()

# End end end
#--------------------------------------------------------------------------------------
