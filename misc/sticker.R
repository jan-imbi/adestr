library(hexSticker)
library(ggplot2)

# Old sticker
imgurl <- here::here("misc/screenshot.png")
sticker <- sticker(imgurl,
  package = "adestr", p_size = 20 * 4000 / 300,
  s_x = 1, s_y = .75, s_width = .9,
  p_color = "black",
  h_fill = "white",
  h_color = "red",
  filename = "man/figures/sticker.png",
  dpi = 4000
)

# New sticker
datadist1 <- Normal(two_armed = FALSE)
H_0 <- PointMassPrior(.0, 1)
H_1 <- PointMassPrior(.4, 1)
ess1 <- ExpectedSampleSize(datadist1, H_1)
power1 <- Power(datadist1, H_1)
toer1  <- Power(datadist1, H_0)
ad_initialD1 <- get_initial_design(
  theta = .4,
  alpha = .025,
  beta  = .2,
  type_design  = "two-stage",
  dist  = datadist1
)
opt1 <- minimize(
  ess1,
  subject_to(
    power1 >= 0.8,
    toer1  <= .025
  ),
  ad_initialD1
)
design1 <- opt1$design
sticker <- sticker(
  plot_sample_mean(Normal(), design1, 0.25, 1, TRUE) +
    scale_x_continuous(name = NULL, breaks=NULL) + scale_y_continuous(name = NULL, breaks =  NULL),
  package = "adestr", p_size = 20 * 400 / 300,
  s_x = 0.9,
  s_y = 0.75,
  s_width = 1.5,
  s_height = 1,
  p_color = "black",
  h_fill = "white",
  h_color = "red",
  filename = "man/figures/sticker.png",
  dpi = 400
)




