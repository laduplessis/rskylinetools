

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

  intervals       <- beastio::getTreeIntervals(tree)
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

skyline_ml <- function(tree, segmentSizes=NULL) {

  if (is.null(segmentSizes)) {
    segmentSizes <- rep(1, tree$Nnode)
  }

  intervals       <- beastio::getTreeIntervals(tree)
  intervals$group <- mapGroups(intervals, segmentSizes, "coalescent")

  popSizes <- rep(0, length(segmentSizes))
  for (i in 1:length(popSizes)) {
      sub   <- (intervals$group == i)
      width <- sum(intervals$length[sub])
      n     <- max(intervals$nlineages[sub])
      k     <- segmentSizes[i]

      popSizes[i] <- width*n*(n-k)/(2*k)
  }
  return(popSizes)
}

#' Epoch Sampling Skyline Plot log-likelihood
#'
#' @export
esp_lk <- function(tree, popSizes, segmentSizes=NULL, sampIntensities, epochSizes=NULL, full=TRUE) {

  if (is.null(segmentSizes)) {
    segmentSizes <- rep(1, tree$Nnode)
  }

  if (is.null(epochSizes)) {
    epochSizes <- length(tree$tip.label)
  }

  intervals         <- beastio::getTreeIntervals(tree)
  intervals$segment <- mapGroups(intervals, segmentSizes, c("coalescent","sample"))
  intervals$epoch   <- mapGroups(intervals, epochSizes, "sample")

  lk <- 0
  for (i in 1:nrow(intervals)) {
    alpha <- choose(intervals$nlineages[i],2)
    beta  <- ifelse(intervals$epoch[i] == 0, 0, sampIntensities[intervals$epoch[i]])
    N     <- popSizes[intervals$segment[i]]

    lk    <- lk - intervals$length[i]*(beta*N + alpha/N)
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
