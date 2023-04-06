library(caracas)
x_1 <- symbol("x_1")
x_2 <- symbol("x_2")
n_1 <- symbol("n_1")
n_2 <- symbol("n_2")
ml_in <- symbol("ml_in")
x2_in <- symbol("x2_in")
mu <- symbol("mu")
mu0 <- symbol("mu0")
mu1 <- symbol("mu1")
n_in <- symbol("n_in")
n2_in <- symbol("n2_in")
c2 <- symbol("c2")
c2_in <- symbol("c2_in")
sigma <- symbol("sigma")
smean2 <- (n_1 * x_1 + n_2 * x_2) / (n_1 + n_2)

left_ml <- (ml_in - mu)
left_lr <- (ml_in - mu) * sqrt(n_in)
left_st <- (ml_in - mu) * n_in
left_sw <-  (x2_in - mu) * sqrt(n2_in) / sigma - c2_in
left_np <- n_in * (2*ml_in*(mu0 - mu1) + mu1*mu1 - mu0*mu0)

right_ml_1 <- (x_1 - mu)
right_ml_2 <- (smean2 - mu)

right_lr_1 <- (x_1 - mu) * sqrt(n_1)
right_lr_2 <- (smean2 - mu) * sqrt(n_1 + n_2)

right_st_1 <- (x_1 - mu) * n_1
right_st_2 <- (smean2 - mu) * (n_1 + n_2)

right_sw <-  (x_2 - mu) * sqrt(n_2) / sigma - c2

right_np_1 <- n_1 * (2*x_1*(mu0 - mu1) + mu1*mu1 - mu0*mu0)
right_np_2 <- (n_1+n_2) * (2*((x_1*n_1 + x_2 * n_2) / (n_1 + n_2)) *(mu0 - mu1) + mu1*mu1 - mu0*mu0)


solve_sys(left_ml, right_ml_1, x_1)
solve_sys(left_ml, right_ml_2, x_2)

solve_sys(left_lr, right_lr_1, x_1)
solve_sys(left_lr, right_lr_2, x_2)

solve_sys(left_st, right_st_1, x_1)
solve_sys(left_st, right_st_2, x_2)

solve_sys(left_st, right_st_1, x_1)
solve_sys(left_st, right_st_2, x_2)

solve_sys(left_sw, right_sw, x_2)

solve_sys(left_np, right_np_1, x_1)
solve_sys(left_np, right_np_2, x_2)








