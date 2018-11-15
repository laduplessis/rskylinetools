# Create test data for Skyline reconstruction tests
rm(list = ls())

#' Function used for generating test data
#'
#' Get sample traces given true time and N vectors, draw <nsamples> piecewise constant
#' approximations of dimension <dimension> with Gaussian noise with standard deviation
#' <noise_sd> added.
getSampleTraces <- function(true_t, true_N, nsamples=10, dimension=5, noise_sd=5, plot=FALSE) {

  sampled_t <- sampled_N <- c()
  for (i in 1:nsamples) {
    sampled   <- c(1,sort(sample(2:(length(true_t)-1),size = dimension-1)),length(true_t))
    sampled_t <- rbind(sampled_t, true_t[sampled])

    mean_N    <- c()
    for (j in 2:length(sampled)) {
      a <- sampled[j-1]
      b <- sampled[j]
      mean_N <- c(mean_N, max(0, round(mean(true_N[a:b])+rnorm(1,0,noise_sd), 0)))
    }

    sampled_N <- rbind(sampled_N, mean_N)
  }

  # Quick plot of sampled skyline traces
  if (plot==TRUE) {
    ssampled_N <- cbind(sampled_N, sampled_N[,ncol(sampled_N)])
    plot(range(sampled_t), range(sampled_N), type='n', ylab='N', xlab='t')
    for (i in 1:nrow(ssampled_N)) {
      lines(sampled_t[i,], ssampled_N[i,], type='s', lty=2)
      #abline(v=sampled_t[i,], lty=3, lwd=0.5)
    }
    grid()
  }

  # Remove column of zeros from times
  return(list(t=sampled_t[,2:ncol(sampled_t)], N=sampled_N))
}



########################################################
# Linear function with only 5 samples and dimension 3
set.seed(25)
maxt   <- 10
true_t <- seq(0,maxt,by=0.01)
true_N <- 20*true_t + 5

linearSample <- getSampleTraces(true_t, true_N, dimension=3, nsamples=5, plot=FALSE)

# On a grid with 5 grid points (calculated by hand)
grid_t     <- seq(0, 10, length.out=5)
linearSampleExpectedGridded <- matrix(c(  41,   9,  59,  27,  20,
                                          41,  49,  59,  27, 124,
                                         110, 136,  59, 140, 124,
                                         176, 136, 171, 140, 124,
                                         176, 136, 171, 140, 124), nrow=5, ncol=5, byrow=FALSE)
colnames(linearSampleExpectedGridded) <- grid_t


# At 5 specific times (calculated by hand)
check_t <- c(1.5, 3, 4.2, 5.4, 7.2)
linearSampleExpectedChecked <- matrix(c(  41,  49,  59,  27,  47,
                                          41, 136,  59,  65, 124,
                                         110, 136,  59, 140, 124,
                                         110, 136,  59, 140, 124,
                                         176, 136, 123, 140, 124), nrow=5, ncol=5, byrow=FALSE)
colnames(linearSampleExpectedChecked) <- check_t


########################################################
# Periodic function with fine grid and a lot of samples
set.seed(25)
maxt   <- 10
true_t <- seq(0,maxt,by=0.01)
true_N <- 50*sin(2*true_t)+100
periodicSample <- getSampleTraces(true_t, true_N, dimension=20, nsamples=1000, plot=FALSE)

grid_t    <- seq(0, 10, length.out=100)
periodicSampleExpectedGridded <- reconstructBackwardSkylineSlow(periodicSample$N, periodicSample$t, grid_t)
colnames(periodicSampleExpectedGridded) <- grid_t


####################################################################
# Exponential function with random time points and a lot of samples
set.seed(25)
maxt   <- 10
true_t <- seq(0,maxt,by=0.01)
true_N <- exp(maxt-0.5*(true_t))
exponentialSample <- getSampleTraces(true_t, true_N, dimension=20, nsamples=1000, plot=FALSE)

randomTimes    <- sort(runif(100, min=0, max=10))
exponentialSampleExpectedRandom <- reconstructBackwardSkylineSlow(exponentialSample$N, exponentialSample$t, randomTimes)
colnames(exponentialSampleExpectedRandom) <- randomTimes

save(linearSample, linearSampleExpectedChecked, linearSampleExpectedGridded, periodicSample, periodicSampleExpectedGridded, exponentialSample, exponentialSampleExpectedRandom, randomTimes,
     file="tests/testDataBackwardSkylines.RData")

