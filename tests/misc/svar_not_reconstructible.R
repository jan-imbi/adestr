#first stage
x1 <- c(1,2)
y1 <- c(1,2)

#second stage
x2 <- c(1,2)
y2 <- c(1,2)

svar1 <- var(x1) + var(y1)
svar2 <- var(x2) + var(y2)
svar <- var(c(x1, x2)) + var(y1, y2)
smean1d <- mean(x1) - mean(y1)
smean2d <- mean(x2) - mean(y2)

c(svar, smean1d, smean2d, svar1, svar2)
###


#first stage
x1 <- c(1,2)
y1 <- c(1,2)

#second stage
x2 <- c(2,3)
y2 <- c(2,3)

svar1 <- var(x1) + var(y1)
svar2 <- var(x2) + var(y2)
svar <- var(c(x1, x2)) + var(y1, y2)
smean1d <- mean(x1) - mean(y1)
smean2d <- mean(x2) - mean(y2)

c(svar, smean1d, smean2d, svar1, svar2)











