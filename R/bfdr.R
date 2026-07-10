bfdr <- function(pips, cutoff) {
  phat <- 1 - pips
  tseq <- seq(min(phat), 1, by=1e-4)
  bfdrt <- numeric(length(tseq))
  counter <- 1
  for(t in tseq) {
    dk <- as.numeric(phat < t)
    bfdrt[counter] <- sum(phat * dk) / sum(dk)
    counter <- counter + 1
  }
  diffs <- abs(cutoff - bfdrt)
  tchoice <- min(which(diffs == min(diffs, na.rm = TRUE)))
  return(1 - tseq[tchoice])
}
