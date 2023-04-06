# This shows why it is a bad idea to integrate over (T1, T_overall) instead of (T1, T2).
q <- qt(.975, 3, 0)
nsim <- 100000
ergmean <- numeric(nsim)
ergvar <- numeric(nsim)
ergmean2 <- list()
ergvar2 <- list()

for (i in 1:nsim){
  a <- rnorm(3)
  t <- t.test(a)
  if (t$statistic > q){
    a <- c(a, rnorm(2))
    ergmean2[[length(ergmean2)+1]] <- mean(a)
    ergvar2[[length(ergvar2)+1]] <- var(a)
  }
  ergmean[[i]] <- mean(a)
  ergvar[[i]] <- var(a)
}
cor.test(ergmean, ergvar)
cor.test(unlist(ergmean2), unlist(ergvar2))
