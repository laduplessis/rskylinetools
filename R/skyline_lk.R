

mapGroups <- function(intervals, groupSizes, eventTypes) {

    events  <- rep(0,nrow(intervals))
    eventNr <- 1
    for (i in 1:nrow(intervals)) {
        events[i] <- eventNr
        if (intervals$nodetype[i] %in% eventTypes) {
            eventNr <- eventNr + 1
        }
    }

    group <- rep(0,nrow(intervals))
    cumGroupSizes <- cumsum(groupSizes)
    for (i in length(cumGroupSizes):1) {
        group[events <= cumGroupSizes[i]] <- i
    }

    return(group)
}


#' Generalized Skyline Plot log-likelihood
#'
#' @export
skyline_lk <- function(tree, popSizes, segmentSizes=NULL, full=TRUE) {

  if (is.null(segmentSizes)) {
      segmentSizes <- rep(1, tree$Nnode)
  }

  intervals       <- getTreeIntervals(tree)
  intervals$group <- mapGroups(intervals, segmentSizes, "coalescent")

  lk <- 0
  for (i in 1:nrow(intervals)) {
      alpha <- choose(intervals$nlineages[i],2)
      N     <- popSizes[intervals$group[i]]

      lk    <- lk - (intervals$length[i]*alpha)/N
      if (intervals$nodetype[i] == "coalescent") {
          lk <- lk - log(N)
          if (full) {
              lk <- lk + log(alpha)
          }
      }
  }

  return(lk)
}

#' Epoch Sampling Skyline Plot log-likelihood
#'
#' @export
esp_lk <- function(tree, segmentSizes, popSizes, epochSizes, sampIntensities) {

  intervals         <- getTreeIntervals(tree)
  intervals$segment <- mapGroups(intervals, segmentSizes, c("coalescent","sample"))
  intervals$epoch   <- mapGroups(intervals, epochSizes, "sample")

  lk <- 0
  for (i in 1:nrow(intervals)) {
    alpha <- choose(intervals$nlineages[i],2)
    beta  <- sampIntensities[intervals$epoch[i]]
    N     <- popSizes[intervals$segment[i]]

    lk    <- lk - intervals$length[i]*(beta*N - alpha/N)
    if (intervals$nodetype[i] == "coalescent") {
      lk <- lk - log(N)
      if (full) {
          lk <- lk + log(alpha)
      }
    } else
    if (intervals$nodetype[i] == "sample") {
      lk <- lk + log(beta*N)
    }
  }

  return(lk)
}


# intervals <- getTreeIntervals(tree)
# groupSizes <- c(2,2,2,2,3)
# groupTypes <- c("sample", "coalescent")
# intervals$group <- mapGroups(intervals, groupSizes, groupTypes)
#
# groupSizes <- c(1,1,1,1,2)
# groupTypes <- c("coalescent")
# intervals$group <- mapGroups(intervals, groupSizes, groupTypes)
#
#
# groupSizes <- c(1,1,1,1,1,1)
# groupTypes <- c("sample")
# intervals$group <- mapGroups(intervals, groupSizes, groupTypes)
#
# tree.test <- read.tree(text = "((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);") # load tree
# sk1 <- skyline(tree.test)
# cat(sprintf("%.12f", sk1$logL))
