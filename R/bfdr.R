bfdr <- function(pips, ctrl) {
    yy <- cumsum(sort(1-pips))/(1:length(pips))
    mm <- yy < ctrl
    idx <- sum(mm)
    if (idx > 0) {
        out <- 1 - sort(1-pips)[idx]
    }
    if (idx == 0) {
        out <- 1
        warning("no valid cutoff found")
    }
    names(out) <- NULL
    return(out)
}
