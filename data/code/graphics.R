# Author: Jan Meis
# Source: https://github.com/jan-imbi/adestr/

# If you did not run the calculations.R script beforehand, you will need to get
# the data from https://github.com/jan-imbi/adestr/tree/master/data
# and save it in the "./data/" directory

# Load all data from the /data/ directory
files <- list.files("data") |>
  grepl(".*\\.rds", x=_) |>
  subset(list.files("data"), subset = _) |>
  (\(x) paste0("data/", x))()
for (f in files){
  assign(stringr::str_extract(f, "(?<=data/).*(?=\\.rds)"), readRDS(f))
}

# All packages and helper functions required to produce the plots
library(adestr)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)
library(Cairo)
sanitize_names <- function(x) {
  if (is.data.frame(x)) {
    xx <- x
    x <- names(xx)
  }
  x <- stringr::str_to_sentence(x)
  x[x=="Mse"] <- "MSE"
  x[x=="Mu"] <- "mu"
  x[x=="Sigma"] <- "sigma"
  if (is.data.frame(xx)) {
    names(xx) <- x
    return(xx)
  } else{
    return(x)
  }
}
sanitize_vars <- function(x) {
  xest <- x$Estimator
  xest[xest=="Bias reduced MLE (iterations=1)"] <- "Bias reduced MLE (1 iteration)"
  xest[xest=="LR test ordering"] <- "LR test ordering CI"
  xest[xest=="MLE ordering"] <- "MLE ordering CI"
  xest[xest=="Score test ordering"] <- "Score test ordering CI"
  xest[xest=="SWCF ordering"] <- "SWCF ordering CI"
  x$Estimator <- xest
  x
}

### MAIN PAPER ###
# Design characteristics: Main paper (1-arm Normal)
f_designs <- adestr:::plot_design(list(designad, designgs), data_distribution = Normal(two_armed = FALSE))
Cairo::CairoTIFF("data/f_designs.tiff", width = 7*2*310, height = 4*2*310, res = 310*2)
f_designs
dev.off()
ggsave("data/f_design.eps" , f_designs, device = "eps", width = 7, height = 4)


# Sample mean: Normal (one-armed)
f_smean <- adestr:::plot_sample_mean(Normal(two_armed = FALSE), designad, mu = seq(-.3, .8, .1), sigma = 1,
                                     combine_components = TRUE, exact = FALSE) + theme(text = element_text(size=9))
Cairo::CairoTIFF("data/f_smean.tiff", width = 7*2*310, height = 4*2*310, res = 310*2)
f_smean
dev.off()
ggsave("data/f_smean.eps" , f_smean, device = "eps", width = 7, height = 4)
saveRDS(f_smean, "data/f_smean.rds")

# Estimator characteristics
tab_est_filter <- tab_est |>
  dplyr::filter(!(estimator%in% c("MedianUnbiasedNeymanPearsonOrdering(mu0=0, mu1=0.4)"))) |>
  sanitize_names() |>
  sanitize_vars()
f_est_mse <- ggplot() +
  geom_line(data = tab_est_filter, aes(x = mu, y=MSE, col = Estimator), linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_bias <- ggplot(data = tab_est_filter, aes(x = mu, y=Bias, col = Estimator)) +
  geom_line(linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est <- ggarrange(f_est_mse, f_est_bias, common.legend = TRUE)
Cairo::CairoTIFF("data/f_est.tiff", width = 7*2*310, height = 4*2*310, res = 310*2)
f_est
dev.off()
ggsave("data/f_est.eps" , f_est, device = "eps", width = 7, height = 4)
saveRDS(f_est, "data/f_est.rds")

# Confidence interval characteristics
tab_ci_filter <- tab_ci |>
  dplyr::filter(!(estimator%in% c(
    "Neyman-Pearson test ordering (mu0=0, mu1=0.4)"
  ))) |>
  sanitize_names() |>
  sanitize_vars()
f_ci_cov <- ggplot(data = tab_ci_filter, aes(x = mu, y=Coverage, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_width <- ggplot(data = tab_ci_filter, aes(x = mu, y=Width, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_test <- ggplot(data = tab_ci_filter, aes(x = mu, y=`Agreement with test decision`, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci <- ggarrange(f_ci_cov, f_ci_width, f_ci_test, common.legend = TRUE, nrow=3)
Cairo::CairoTIFF("data/f_ci.tiff", width = 5.6*2*310, height = 7.6*2*310, res = 310*2)
f_ci
dev.off()
ggsave("data/f_ci.eps" , f_ci, device = "eps", width = 5.6, height = 7.6)
saveRDS(f_ci, "data/f_ci.rds")

# P-values / stage-wise ordering representation for an adaptive design
f_p_sw <- adestr:::plot_p(StagewiseCombinationFunctionOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[1], subdivisions = 200L)
f_p_ml <- adestr:::plot_p(MLEOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[2], subdivisions = 200L)
f_p_lr <- adestr:::plot_p(LikelihoodRatioOrderingPValue(), Normal(two_armed = FALSE),  designad, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)
f_p_st <- adestr:::plot_p(ScoreTestOrderingPValue(), Normal(two_armed = FALSE), designad, 0, 1, boundary_color = scales::hue_pal()(5)[4], subdivisions = 200L)
f_p_np <- adestr:::plot_p(NeymanPearsonOrderingPValue(), Normal(two_armed = FALSE), designad, -1, 1, boundary_color = scales::hue_pal()(5)[5], subdivisions = 200L)
f_p_reject <- adestr:::plot_rejection_regions(
  list(StagewiseCombinationFunctionOrderingPValue(), MLEOrderingPValue(),  LikelihoodRatioOrderingPValue(), ScoreTestOrderingPValue(), NeymanPearsonOrderingPValue()),
  Normal(two_armed = FALSE), designad, 0, 1, subdivisions = 200L)
f_p <- ggarrange(f_p_sw, f_p_ml, f_p_lr, f_p_st, f_p_np, f_p_reject,  ncol = 2, nrow = 3,
                 legend.grob = gridExtra::gtable_rbind(get_legend(f_p_reject), get_legend(f_p_sw)))
Cairo::CairoTIFF("data/f_p.tiff", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
ggsave("data/f_p.eps", f_p, width = 900*3.5/310, height = 700*3.5/310)
dev.off()
Cairo::CairoPNG("data/f_p.png", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
dev.off()

### SUPPLEMENT ###
# P-values / stage-wise ordering representation for a group-sequential design
f_p_sw <- adestr:::plot_p(StagewiseCombinationFunctionOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[1], subdivisions = 200L)
f_p_ml <- adestr:::plot_p(MLEOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[2], subdivisions = 200L)
f_p_lr <- adestr:::plot_p(LikelihoodRatioOrderingPValue(), Normal(two_armed = FALSE),  designgs, 0, 1, boundary_color = scales::hue_pal()(5)[3], subdivisions = 200L)
f_p_st <- adestr:::plot_p(ScoreTestOrderingPValue(), Normal(two_armed = FALSE), designgs, 0, 1, boundary_color = scales::hue_pal()(5)[4], subdivisions = 200L)
f_p_np <- adestr:::plot_p(NeymanPearsonOrderingPValue(), Normal(two_armed = FALSE), designgs, -1, 1, boundary_color = scales::hue_pal()(5)[5], subdivisions = 200L)
f_p_reject <- adestr:::plot_rejection_regions(
  list(StagewiseCombinationFunctionOrderingPValue(), MLEOrderingPValue(),  LikelihoodRatioOrderingPValue(), ScoreTestOrderingPValue(), NeymanPearsonOrderingPValue()),
  Normal(two_armed = FALSE), designgs, 0, 1, subdivisions = 200L)
f_p <- ggarrange(f_p_sw, f_p_ml, f_p_lr, f_p_st, f_p_np, f_p_reject,  ncol = 2, nrow = 3,
                 legend.grob = gridExtra::gtable_rbind(get_legend(f_p_reject), get_legend(f_p_sw)))
Cairo::CairoTIFF("data/f_p_gs.tiff", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
dev.off()
Cairo::CairoPNG("data/f_p_gs.png", width = 900*3.5*2, height = 700*3.5*2, res = 310*2)
f_p
dev.off()

# Design characteristics: Appendix (2-arm t-distribution)
f_designs_t <-  adestr:::plot_design(list(designad_t, designgs_t), data_distribution = Student(two_armed = TRUE))
f_designs_t
saveRDS(f_designs_t, "data/f_designs_t.rds")

# Sample mean of treatment group: T (two-armed)
f_smean_t <- adestr:::plot_sample_mean(Student(two_armed = TRUE), designad_t, mu = seq(-.3, .8, .1), sigma = 1,
                                       combine_components = TRUE, exact = FALSE)
f_smean_t
saveRDS(f_smean_t, "data/f_smean_t.rds")

# Estimator characteristics
tab_est_filter_t <- tab_est_t |>
  dplyr::filter(!(estimator%in% c("MedianUnbiasedNeymanPearsonOrdering(mu0=0, mu1=0.4)"))) |>
  sanitize_names() |>
  sanitize_vars()
f_est_mse_t <- ggplot() +
  geom_line(data = tab_est_filter_t, aes(x = mu, y=MSE, col = Estimator), linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_bias_t <- ggplot(data = tab_est_filter_t, aes(x = mu, y=Bias, col = Estimator)) +
  geom_line(linewidth = .5) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr() +
  theme(text = element_text(size=9)) +
  guides(color=guide_legend(ncol=3))
f_est_t <- ggarrange(f_est_mse_t, f_est_bias_t, common.legend = TRUE)
f_est_t
saveRDS(f_est_t, "data/f_est_t.rds")

# Confidence interval characteristics
tab_ci_filter_t <- tab_ci_t |>
  dplyr::filter(!(estimator%in% c(
    "Neyman-Pearson test ordering (mu0=0, mu1=0.4)"
  ))) |>
  sanitize_names() |>
  dplyr::filter(-0.74 < mu & mu < 1.33)
f_ci_cov_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=Coverage, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_width_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=Width, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_test_t <- ggplot(data = tab_ci_filter_t, aes(x = mu, y=`Agreement with test decision`, col = Estimator)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(Design)) +
  scale_x_continuous(TeX("$\\mu$")) +
  theme_pubr()
f_ci_t <- ggarrange(f_ci_cov_t, f_ci_width_t, f_ci_test_t, common.legend = TRUE, nrow=3)
f_ci_t
saveRDS(f_ci_t, "data/f_ci_t.rds")
