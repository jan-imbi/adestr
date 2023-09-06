# Helper function to plot designs
#' @importFrom scales percent
plot_design <- function(design, data_distribution = Normal(two_armed = FALSE)){
  design <- TwoStageDesignWithCache(design)
  two_armed <- data_distribution@two_armed
  if (is(data_distribution, "Student")) {
    z1tex <- TeX("$t_1$")
    c2tex <- TeX("$c_2(t_1)$")
  } else {
    z1tex <- TeX("$z_1$")
    c2tex <- TeX("$c_2(z_1)$")
  }
  contl <- max(sapply(design, \(x) x@c1e - x@c1f))
  minc1f <- min(sapply(design, \(x) x@c1f))
  maxc1e <- max(sapply(design, \(x) x@c1e))
  const <- (1-2/(1+sqrt(5)))
  plotdata <- list()
  n1_dodge <- seq(-.1, .1, length.out = length(design))
  for (i in seq_along(design)) {
    d <- design[[i]]
    x1 <- seq(minc1f - const * contl, d@c1f, length.out = 200)
    x2 <- seq(d@c1f, d@c1e, length.out = 200)
    x3 <- seq(d@c1e, maxc1e + const * contl, length.out = 200)
    n1 <- n1(d, round=FALSE) + n1_dodge[[i]]
    n2 <- n2_extrapol(d, x2)
    c2 <- c2_extrapol(d, x2)
    cp <- pnorm(c2_extrapol(d, x2), mean = 0.4*sqrt(n2_extrapol(d, x2) / (1L + two_armed)), lower.tail = FALSE)
    plotdata[[length(plotdata) + 1L]] <- data.frame(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      n1 = n1,
      n2 = n2,
      c2 = c2,
      cp = cp,
      label = toString(d)
    )
  }
  dat <- do.call("rbind", plotdata)
  pltn <- pltc2 <- pltcp <-  ggplot(data = dat)
  pltn <- pltn +
    geom_line(aes(x = x1, y = n1,      color = .data$`label`), linewidth=1) +
    geom_line(aes(x = x2, y = n1 + n2, color = .data$`label`), linewidth=1) +
    geom_line(aes(x = x3, y = n1,      color = .data$`label`), linewidth=1)
  pltc2 <- pltc2 +
    geom_line(aes(x = x2, y = c2, color = .data$`label`), linewidth=1) +
    geom_line(aes(x = x2, y = c2, color = .data$`label`), linewidth=1)
  pltcp <- pltcp +
    geom_line(aes(x = x2, y = cp, color = .data$`label`), linewidth=1) +
    geom_line(aes(x = x2, y = cp, color = .data$`label`), linewidth=1)
  pltn <- pltn +
    scale_x_continuous(name = z1tex, breaks = unique(round(x2)) ) +
    scale_y_continuous(name = "Overall sample size", limits = c(0, 10*ceiling(max((dat$n2+dat$n1) /10)) + 10 * (1L + two_armed)) ) +
    theme_pubr() +
    theme(text = element_text(size=15)) +
    labs(color = "Type of design")
  pltc2 <- pltc2 +
    scale_x_continuous(name = z1tex, breaks = unique(round(x2)))+
    scale_y_continuous(name = c2tex) +
    theme_pubr() +
    theme(text = element_text(size=15)) +
    labs(color = "Type of design")
  pltcp <- pltcp +
    scale_x_continuous(name = z1tex, breaks = unique(round(x2)))+
    scale_y_continuous(name = "Conditional power",
                       labels = scales::percent) +
    theme_pubr() +
    theme(text = element_text(size=15)) +
    labs(color = "Type of design")
  ggarrange(pltn, pltc2, pltcp, ncol=3, common.legend = TRUE)
}
