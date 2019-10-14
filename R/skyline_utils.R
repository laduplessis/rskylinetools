

#' Function used for testing
#'
#' Get sample traces given true time and N vectors, draw <nsamples> piecewise constant
#' approximations of dimension <dimension> with Gaussian noise with standard deviation
#' <noise_sd> added.
#'
#' @export
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




#'     Case 1:
#'         skyline   = m x n matrix
#'         times     = m x n matrix
#'         gridtimes = n-dim vector
#'
#'         Case 1a:
#'             times[i,1] < gridtimes[1] for all i     -> Add column of Inf to the right of times
#'         Case 1b:
#'             times[i,1] >= gridtimes[1] for some i
#'             times[i,n] >  gridtimes[n] for all i    -> Add column of -Inf to the left of times
#'         Case 1c:
#'             times[i,1] >= gridtimes[1] for some i   -> Add column of -Inf to the left of times
#'             times[i,n] <= gridtimes[n] for some i   -> Change last column of times to Inf
#'
#'     Case 2:
#'         skyline   = m x n matrix
#'         times     = m x n-1 matrix
#'         gridtimes = n-dim vector
#'
#'         times[i,1] > gridtimes[1] for all i
#'         times[i,n-1] < gridtimes[n] for all i      -> Add column of -Inf to the left and Inf to the right of times
#'
#'     Anything else:
#'         Error!


#' @export
getChangeTimes <- function(trees, groupsizes, eventtypes="coalescent") {

  if (nrow(groupsizes) != length(trees)) stop("Error! Different number of trees to number of states in log file")

  # Calculate change-point times from logged trees and logged group sizes
  changetimes <- matrix(0,ncol=ncol(groupsizes),nrow=nrow(groupsizes))
  for (i in 1:nrow(groupsizes)) {
    events          <- beastio::getTreeIntervals(trees[[i]])
    coaltimes       <- events$height[events$nodetype == "coalescent"]
    samplingtimes   <- events$height[events$nodetype == "sample"]
    eventtimes      <- events$height

    groupshifts     <- cumsum(simplify2array(groupsizes[i,]))

    if (eventtypes == "coalescent") {
        changetimes[i,] <- coaltimes[groupshifts]
    } else
    if (eventtypes == "sampling") {
        changetimes[i,] <- samplingtimes[groupshifts]
    } else {
        changetimes[i,] <- eventtimes[groupshifts]
    }
  }

  return(changetimes)
}




# Rewrite with findInterval...

#' Reconstructs a skyline that is measured backwards-in-time (t = 0 is the present)
#' THIS IS THE SLOW AND CAREFUL FUNCTION THAT ONLY USES PRIMITIVES AND NO VECTORIZATION!
#'
#' For state i and interval j, skyline[i,j] is the estimate between times[i,j-1] and times[i,j].
#' Thus,
#'
#'     skyline   = m x n matrix
#'     times     = m x n matrix
#'     gridtimes = n-dim vector
#'
#'     times[i] < times[i+1] (sorted)
#'     times[i-1] < times[i] (sorted)
#'     times[i,j] != 0
#'
#'     skyline[i,1] <=> t < times[i,1]                  j == 1 (implicitly set t_0 to 0)
#'     skyline[i,j] <=> times[i,j-1] < t <= times[i,j]  j > 1
#'     skyline[i,j] <=> t > times[i,j]                  j == ncol(skyline)
#'
#'
#' @param skyline    : Skyline estimates in a matrix or data.frame form (rows = MCMC states, cols = skyline estimates during intervals)
#' @param times      : Times associated with the end-times of each skyline interval (rows = MCMC states, cols = interval start-times)
#'                     Must be the same dimension as skyline.
#' @param maxtime    : Maximum time of the skyline function. If this is a scalar it is assumed that the same maxtime is used for every row in skyline
#' @param gridtimes  : The regular times to reconstruct the marginal posterior on.
#'
#' @export
reconstructBackwardSkylineSlow <- function(skyline, times, gridtimes) {

  n <- ncol(skyline)   # skyline dimension
  m <- nrow(skyline)   # number of samples (states)

  # Check skyline times dimensions
  if (nrow(times) != m || ncol(times) != n) stop("Error! Wrong number of rows or columns in skyline times matrix")

  # Check skyline times order
  for (i in 1:m) {
    if (is.unsorted(times[i,])) stop(paste0("Error! Skyline times row ",i," is not sorted"))
  }

  # Check skyline times first column = 0
  # if (sum(times[,1] != 0) > 0) stop("Error! First column of skyline times must be equal to 0")

  # Check that all times are > 0
  if (min(times) <= 0) stop("Error! All times must be > 0")

  # Set first and last times to -Inf and Inf (so all times in gridtimes fit in)
  times[,n] <- rep(Inf, m)
  times     <- cbind(rep(-Inf,m), times)

  # Ensure gridtimes are sorted
  gridsize  <- length(gridtimes)
  gridtimes <- sort(gridtimes)
  gridrange <- range(gridtimes)

  skyline_gridded <- matrix(0, nrow=m, ncol=gridsize)
  for (i in 1:m) {

    j <- 1  # Only need to reset index at start of each sample
    for (k in 1:gridsize) {
      # Find interval that gridtimes[k] falls in
      while (times[i,j+1] < gridtimes[k]) {
        j <- j + 1
      }

      # skyline[i,j] <=> times[i,j] < t <= times[i,j+1]  j < ncol(skyline)
      # skyline[i,j] <=> times[i,j-1] < t <= times[i,j]  j > 1

      if (times[i,j] < gridtimes[k] && times[i,j+1] >= gridtimes[k]) {
        skyline_gridded[i,k] <- skyline[i,j]
      } else {
        stop(paste0("Error! Problem at sample ",i,", interval ",j," grid cell ",k,": ", gridtimes[k], " not in (",times[i,j],",",times[i,j+1],"]"))
      }
    }
  }

  colnames(skyline_gridded) <- gridtimes
  return(skyline_gridded)
}

