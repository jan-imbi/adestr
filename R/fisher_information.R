#nocov start
fisher_information <- function(design, mu, sigma, two_armed) {
  sapply(mu, \(x)
    integrate_over_sample_space(
      Normal(two_armed = two_armed),
      use_full_twoarm_sampling_distribution = FALSE,
      design = design,
      g1 = \(smean1, n1, mu, sigma, two_armed, ...) (n1 / (sigma * (1L + two_armed)) * (mu - smean1))^2 ,
      g2 = \(smean1, smean2, n1, n2, mu, sigma, two_armed, ...) ((n1 + n2) / (sigma * (1L + two_armed)) * (mu - (smean1*n1 + smean2 * n2) / (n1 + n2)))^2,
      mu = x,
      sigma = sigma
      )$overall_integral$integral
  )
}

make_inverse_fisher_information_table <- function(designs, mu, sigma, two_armed) {
  res_list <- list()
  for (i in seq_along(designs)){
    design <- designs[[i]]
    InverseFisherInformation <- 1/fisher_information(design, mu, sigma, two_armed)
    res_list[[i]] <- data.frame(Design= attr(design, "label"),
                                InverseFisherInformation = InverseFisherInformation,
                                mu = mu)
  }
  do.call(rbind, res_list)
}
#nocov end
