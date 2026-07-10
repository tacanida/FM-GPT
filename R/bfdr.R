bfdr <- function(pips, alpha = 0.1) {
    
    phat <- sort(1 - pips)
    
    bfdr <- cumsum(phat) / seq_along(phat)
    
    k <- max(c(0, which(bfdr <= alpha)))
    
    if(k == 0)
        return(1)
    
    1 - phat[k]
}
