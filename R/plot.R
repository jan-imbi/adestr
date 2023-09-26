#' Plot performance scores for point and interval estimators
#'
#' This function extract the values of mu and the score values and a facet plot with
#' one facet per score. If the input argument is a list, the different estimators
#' will be displayed in the same facets, differentiated by color.
#'
#' @param x an output object from evaluate_estimator (\code{EstimatorScoreResult}) or a list
#' of such objects (\code{EstimatorScoreResultList}).
#' @param y unused.
#' @param ... additional arguments handed down to ggplot.
#' @export
#' @importFrom ggplot2 ggplot scale_x_continuous geom_line facet_wrap
#' @importFrom latex2exp TeX
#' @returns a \code{\link{ggplot2}} object visualizing the score values.
#' @examples
#' score_result1 <- evaluate_estimator(
#'   MSE(),
#'   estimator = SampleMean(),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu=seq(-.75, 1.32, 0.03),
#'   sigma=1)
#' # Plotting the result of evaluate_estimator
#' plot(score_result1)
#'
#' score_result2 <- evaluate_estimator(
#'   MSE(),
#'   estimator = AdaptivelyWeightedSampleMean(w1 = 0.8),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu=seq(-.75, 1.32, 0.03),
#'   sigma=1)
#' # Plotting a list of different score results
#' plot(c(score_result1, score_result2))
setMethod("plot", signature = "EstimatorScoreResultList", definition =
            function(x, ...) {
              dat <- data.frame()
              for (estimatorscoreresult in x) {
                plot_list <- names(estimatorscoreresult@results)
                for (score in plot_list){
                  dat <- rbind(dat,
                               data.frame(mu = estimatorscoreresult@mu,
                                          Score = estimatorscoreresult@results[[score]],
                                          Estimator = toString(estimatorscoreresult@estimator),
                                          score_name = score
                               )
                  )
                }
              }
              ggplot(data = dat, mapping = aes(x = .data$mu, y = .data$Score, col = .data$Estimator), ...) +
                scale_x_continuous(name = TeX("$\\mu$")) +
                geom_line() +
                facet_wrap(vars(.data$score_name), scales = "free_y")
            })
#' @inherit plot,EstimatorScoreResultList-method
setMethod("plot", signature = "EstimatorScoreResult", definition =
            function(x, ...) {
              l <- EstimatorScoreResultList(x)
              plot(l, ...)
            })
#' @inherit plot,EstimatorScoreResultList-method
#' @importFrom graphics plot.default
setMethod("plot", signature = "list", definition =
            function(x, ...) {
              if (is(x[[1]], "EstimatorScoreResult")) {
                class(x) <- c("EstimatorScoreResultList", class(x))
                plot(x, ...)
              } else {
                plot.default(x, ...)
              }
            })

#' Plot p-values and implied rejection boundaries
#'
#' Creates a plot of the p-values and implied rejection boundaries on a grid
#' of values for the first and second-stage test statistics.
#'
#' When the first-stage test statistic lies below the futility threshold (c1f) or
#' above the early efficacy threshold (c1e) of the \code{TwoStageDesign},
#' there is no second-stage test statistics. The p-values in these regions are only
#' based on the first-stage values.
#' For first-stage test statistic values between c1f and c1e, the first and second-stage
#' test statistic determine the p-value.
#'
#' The rejection boundary signals the line where
#'
#' @inheritParams evaluate_estimator
#' @param boundary_color color of the implied rejection boundary.
#' @param subdivisions number of subdivisions per axis for the grid of test statistic values.
#' @param ... additional arguments handed down to ggplot
#'
#' @returns a \code{\link{ggplot2}} object visualizing the p-values on a grid of possible test-statistic values.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_tile geom_line geom_segment scale_color_manual scale_fill_gradient scale_x_continuous
#' @importFrom ggpubr theme_pubclean
#' @importFrom latex2exp TeX
#'
#' @examples
#' plot_p(estimator = StagewiseCombinationFunctionOrderingPValue(),
#'   data_distribution = Normal(FALSE),
#'   design = get_example_design(),
#'   mu = 0,
#'   sigma = 1)
plot_p <- function(estimator, data_distribution, design, mu = 0, sigma, boundary_color="lightgreen", subdivisions = 100, ...){
  design <- TwoStageDesignWithCache(design)
  two_armed <- data_distribution@two_armed
  n1 <- n1(design, round = FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  contl <- design@c1e - design@c1f
  minx <- design@c1f - 1.8*(1-2/(1+sqrt(5))) * contl
  maxx <- design@c1e + 2.2*(1-2/(1+sqrt(5))) * contl
  miny <- design@c1f - 4*(1-2/(1+sqrt(5))) * contl
  maxy <- design@c1e + (1-2/(1+sqrt(5))) * contl

  gridx <- seq(minx, maxx, length.out = subdivisions)
  gridy <- seq(miny, maxy, length.out = subdivisions)
  region <- expand.grid(T2 = gridy,
                        T1 = gridx)

  stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                   data_distribution =  data_distribution,
                                                   use_full_twoarm_sampling_distribution = FALSE,
                                                   design = design, sigma = sigma, exact = FALSE)
  p1 <- stagewise_estimators[[1L]]
  p2 <- stagewise_estimators[[2L]]

  p_futility <- sapply(gridx[gridx < design@c1f], \(T1) p1(design = design, smean1 = T1 * se1, n1 = n1, sigma = sigma, two_armed = two_armed))
  p_efficacy<- sapply(gridx[gridx > design@c1e], \(T1) p1(design = design, smean1 = T1 * se1, n1 = n1, sigma = sigma, two_armed = two_armed))

  continuation_region <- region[design@c1f <= region$T1 & region$T1 <= design@c1e,]
  p_continuation <- mapply(\(T1, T2) {
    n2 <- n2_extrapol(design, T1)
    se2 <- sigma_to_se(sigma, n2, two_armed)
    ret <- p2(
      design = design,
      smean1 = T1 * se1,
      smean2 = T2 * se2,
      n1 = n1,
      n2 = n2,
      sigma = sigma,
      two_armed = two_armed
    )
  },
  T1 = continuation_region$T1,
  T2 = continuation_region$T2)

  continuation_x_c2 <- seq(design@c1f, design@c1e, .01)
  alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)
  if (is(estimator, "StagewiseCombinationFunctionOrderingPValue")){
    continuation_y_c2 <- c2_extrapol(design, continuation_x_c2)
  } else{
    continuation_y_c2 <- implied_c2(design, continuation_x_c2, p2, sigma, two_armed, alpha)
  }

  draw_line_3 <- FALSE
  if (p1(design = design, smean1 = z_to_smean(design@c1f, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) > alpha) {
    yend1 <- maxy
  } else{
    yend1 <- miny
    draw_line_3 <- TRUE
    x_line_3 <- uniroot(\(x){
      p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
    },
    lower = -4, upper = 4, extendInt="yes")$root
  }
  if (p1(design = design, smean1 = z_to_smean(design@c1e, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) < alpha) {
    yend2 <- miny
  } else{
    yend2 <- maxy
    draw_line_3 <- TRUE
    x_line_3<- uniroot(\(x){
      p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
    },
    lower = -4, upper = 4, extendInt="yes")$root
  }
  continuation_x_c2[c(1, length(continuation_x_c2))] <- continuation_x_c2[c(1, length(continuation_x_c2))]
  c2_comb <- data.frame(x = c(continuation_x_c2[1L],  continuation_x_c2, continuation_x_c2[length(continuation_x_c2)]),
                        y = c(yend1, continuation_y_c2, yend2),
                        PValue = estimator@label)
  p_comb <- rbind(
    cbind(
      T1 = gridx[gridx < design@c1f],
      T2 = rep(0, length(gridx[gridx < design@c1f])),
      p = p_futility
    ),
    cbind(continuation_region, p = p_continuation),
    cbind(
      T1 = gridx[gridx > design@c1e],
      T2 = rep(0, length(gridx[gridx > design@c1e])),
      p = p_efficacy
    )
  )

  p_comb <- cbind(region,
                  p = c(rep(p_futility, each = length(gridx)),
                        p_continuation,
                        rep(p_efficacy, each = length(gridx))
                  ))
  limitsx <- c(min(gridx - .3), max(gridx + .3))
  limitsy <- c(min(gridy - .3), max(gridy + .3))
  if (is(data_distribution, "Normal")) {
    xtxt <- TeX("$z_1$")
    ytxt <- TeX("$z_2$")
  } else if (is(data_distribution, "Student")) {
    xtxt <- TeX("$t_1$")
    ytxt <- TeX("$t_2$")
  } else {
    stop("unsupported data distribution")
  }
  plt <- ggplot(...) +
    geom_tile(data = p_comb, aes(x = .data$`T1`, y = .data$`T2`, fill = .data$`p`)) +
    geom_line(data = c2_comb, aes(x = .data$x, y = .data$y), color = boundary_color, linewidth=1) +
    scale_color_manual() +
    scale_fill_gradient(low="blue", high="red") +
    scale_x_continuous(name =xtxt,
                       limits = limitsx,
                       breaks = unique(round(seq(limitsx[1], limitsx[2], .5) )) ) +
    scale_y_continuous(name =ytxt,
                       limits = limitsy,
                       breaks = unique(round(seq(limitsy[1], limitsy[2], .5) ))) +
    theme_pubclean() +
    theme(plot.title = element_text(hjust = 0.5), legend.key.width = unit(1, 'cm')) +
    ggtitle(TeX(toTeX(estimator)))
  if (draw_line_3) {
    plt <- plt + geom_segment(data = c2_comb, mapping = aes(x=x_line_3, xend=x_line_3, y=miny, yend = maxy), color = boundary_color, linewidth=1)
  }
  plt
}

### Some other (unexported) plotting methods used in the paper ###

#' @importFrom latex2exp TeX
plot_rejection_regions <- function(estimators, data_distribution, design, mu, sigma,  subdivisions = 100, ...){
  design <- TwoStageDesignWithCache(design)
  two_armed <- data_distribution@two_armed
  n1 <- n1(design, round = FALSE)
  se1 <- sigma_to_se(sigma, n1, two_armed)
  contl <- design@c1e - design@c1f
  minx <- design@c1f - 1.8*(1-2/(1+sqrt(5))) * contl
  maxx <- design@c1e + 2.2*(1-2/(1+sqrt(5))) * contl
  miny <- design@c1f - 4*(1-2/(1+sqrt(5))) * contl
  maxy <- design@c1e + (1-2/(1+sqrt(5))) * contl

  gridx <- seq(minx, maxx, length.out = subdivisions)
  gridy <- seq(miny, maxy, length.out = subdivisions)
  region <- expand.grid(T2 = gridy,
                        T1 = gridx)
  continuation_region <- region[design@c1f <= region$T1 & region$T1 <= design@c1e,]

  continuation_x_c2 <- seq(design@c1f, design@c1e, .01)
  alpha <- adoptr_alpha_shifted_design_kv(design, 0, 0, 0)

  c2_comb_list <- list()
  draw_line_3 <- FALSE
  for (i in seq_along(estimators)) {
    estimator <- estimators[[i]]
    stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                     data_distribution =  data_distribution,
                                                     use_full_twoarm_sampling_distribution = FALSE,
                                                     design = design, sigma = sigma, exact = FALSE)
    p1 <- stagewise_estimators[[1L]]
    p2 <- stagewise_estimators[[2L]]
    if (is(estimator, "StagewiseCombinationFunctionOrderingPValue")){
      continuation_y_c2 <- c2_extrapol(design, continuation_x_c2)
    } else{
      continuation_y_c2 <- implied_c2(design, continuation_x_c2, p2, sigma, two_armed, alpha)
    }
    x_line_3 <- NA
    if (p1(design = design, smean1 = z_to_smean(design@c1f, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) > alpha) {
      yend1 <- maxy
    } else{
      yend1 <- miny
      draw_line_3 <- TRUE
      x_line_3 <- uniroot(\(x){
        p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
      },
      lower = -4, upper = 4, extendInt="yes")$root
    }
    if (p1(design = design, smean1 = z_to_smean(design@c1e, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) < alpha) {
      yend2 <- miny
    } else{
      yend2 <- maxy
      draw_line_3 <- TRUE
      x_line_3<- uniroot(\(x){
        p1(design = design, smean1 = z_to_smean(x, n1, sigma, two_armed), n1 = n1, sigma = sigma, two_armed = two_armed) - alpha
      },
      lower = -4, upper = 4, extendInt="yes")$root
    }
    continuation_x_c2[c(1, length(continuation_x_c2))] <- continuation_x_c2[c(1, length(continuation_x_c2))]
    c2_comb <- data.frame(x = c(continuation_x_c2[1L],  continuation_x_c2, continuation_x_c2[length(continuation_x_c2)]),
                          y = c(yend1, continuation_y_c2, yend2),
                          x_line_3 = x_line_3,
                          Ordering = if (is(estimator, "NeymanPearsonOrderingPValue")) "Neyman-Pearson test ordering" else  substr(toTeX(estimator), 1, nchar(toTeX(estimator)) -8))
    c2_comb_list[[i]] <- c2_comb

  }
  c2_comb <- do.call("rbind", c2_comb_list)
  c2_comb$Ordering <- factor(c2_comb$Ordering, levels = unique(c2_comb$Ordering))

  limitsx <- c(min(gridx - .3), max(gridx + .3))
  limitsy <- c(min(gridy - .3), max(gridy + .3))
  plt <- ggplot() +
    geom_line(data = c2_comb, aes(x = .data$x, y = .data$y, color = .data$Ordering), linewidth=1) +
    scale_color_discrete(labels = sapply(as.character(unique(c2_comb$Ordering)), TeX)) +
    scale_x_continuous(name =TeX("$z_1$"),
                       limits = limitsx,
                       breaks = unique(round(seq(limitsx[1], limitsx[2], .5) )) ) +
    scale_y_continuous(name =TeX("$z_2$"),
                       limits = limitsy,
                       breaks = unique(round(seq(limitsy[1], limitsy[2], .5) ))) +
    theme_pubclean() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All rejection boundaries combined")
  if (draw_line_3) {
    c2_comb2 <- c2_comb[!is.na(c2_comb$x_line_3),]
    plt <- plt + geom_segment(data = c2_comb2, mapping = aes(x=x_line_3, xend=x_line_3, y=miny, yend = maxy, color=.data$Ordering), linewidth=1)
  }
  plt
}

setGeneric("plot_sample_mean", \(data_distribution, design, mu, sigma, combine_components = TRUE, plot_treatment_group_if_twoarm = FALSE,
                                 p_limit = .0001, subdivisions = 100L,
                                 exact = FALSE, facets = 3L, ...) standardGeneric("plot_sample_mean"))
#' @import ggplot2
#' @importFrom forcats as_factor
setMethod("plot_sample_mean", signature("DataDistribution", "TwoStageDesign"),
          \(data_distribution, design, mu, sigma, combine_components, p_limit, subdivisions, exact, facets, ...) {
            two_armed <- data_distribution@two_armed
            n1 <- n1(design, round = FALSE)
            se1 <- sigma_to_se(sigma, n1, two_armed)
            smean_list <- list()
            for (m in mu){
              if (is(data_distribution, "Student")) {
                gridx <- suppressWarnings(seq(qt(1-p_limit, df = n1, ncp = m/se1, lower.tail = FALSE),
                                              qt(p_limit, df = n1, ncp = m/se1, lower.tail = FALSE),
                                              length.out = subdivisions)) * se1
                gridx <- gridx[!(abs(gridx)<.02)]
              } else {
                gridx <- seq(qnorm(1-p_limit, mean = m, sd = se1, lower.tail = FALSE),
                             qnorm(p_limit, mean = m, sd = se1, lower.tail = FALSE),
                             length.out = subdivisions)
              }
              if (plot_treatment_group_if_twoarm){
                y <- dsmeanT(data_distribution, design, smeanT = gridx, mu = m, sigma = sigma,
                             combine_components = combine_components, exact = exact)
              } else{
                y <- dsmean(data_distribution, design, smean = gridx, mu = m, sigma = sigma, two_armed = two_armed,
                            combine_components = combine_components, exact = exact)
              }

              if (combine_components)
                smean_list[[length(smean_list)+1L]] <- data.frame(smean = gridx, Density = y, mu = (paste0("mu == ",format(round(m, digits = 2)))))
            }
            if (combine_components){
              smean_dat <- do.call("rbind", smean_list)
              smean_dat$mu <- as_factor(smean_dat$mu)
              plt <- ggplot(data = smean_dat, aes(x = .data$`smean`, y = .data$`Density`)) +
                geom_line(size = 1) +
                scale_x_continuous(name = "Overall sample mean") +
                theme_pubclean()
              if (length(mu) > 1L) {
                plt <- plt + facet_wrap(vars(mu), labeller = label_parsed)
              }
              return(plt)
            } else {
              if (exact) {
                y[["0"]] <- y[["futility"]] + y[["efficacy"]]
                y$efficacy <- NULL
                y$futility <- NULL
                ys <- c(y["0"], y[seq(1L, length(y)-1L, length.out = facets)])
              }
              else
                ys <- c(y["futility"], y["continuation"], y["efficacy"])
              smean_dat <- data.frame(smean = rep(gridx, times = length(ys)),
                                      Density = do.call(c, ys),
                                      n = if (exact) paste0("n2 = ",rep(names(ys), each = length(gridx))) else rep(c("futility", "continuation", "efficacy"), each = length(gridx)) )
              smean_dat$n <- as_factor(smean_dat$n)
              plt <- ggplot(data = smean_dat, aes(x = .data$`smean`, y = .data$`Density`)) +
                geom_line(size = 1) +
                scale_x_continuous(name = "Sample mean") +
                facet_wrap(vars(.data$n))
              return(plt)
            }
          })

# Helper function to plot designs
#' @importFrom scales percent
#' @importFrom ggpubr theme_pubr ggarrange
#' @import ggplot2
plot_design <- function(design, data_distribution = Normal(two_armed = FALSE)){
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
    d <- TwoStageDesignWithCache(d)
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
