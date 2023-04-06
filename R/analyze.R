Results <- setClass("Results", slots = c(data ="data.frame",
                                         summary_data = "data.frame",
                                         data_distribution = "DataDistribution",
                                         design = "TwoStageDesign",
                                         sigma = "numeric",
                                         results = "list"))
#' Analyze a dataset
#'
#' @param data a data.frame containing the data to be analyzed
#' @inheritParams evaluate_estimator
#'
#' @return \code{Results} object containing the values of the estimators
#' when applied to data.
#' @export
#'
#' @examples
#' set.seed(321)
#' dat <- data.frame(
#'   endpoint = c(rnorm(28, .2, 1), rnorm(28, 0, 1),
#'                rnorm(23, .2, 1), rnorm(23, 0, 1)),
#'   group = factor(rep(c("ctl", "trt", "ctl", "trt"),
#'                      c(28,28,23,23))),
#'   stage = rep(c(1L, 2L), c(56, 46))
#' )
#' analyze(
#'   data = dat,
#'   estimator = c(get_example_estimators()),
#'   data_distribution = Normal(),
#'   sigma = 1,
#'   design = get_example_design()
#' )
setGeneric("analyze", function(data,
                               estimator,
                               data_distribution,
                               use_full_twoarm_sampling_distribution = FALSE,
                               design,
                               sigma,
                               exact = FALSE) standardGeneric("analyze"))
#' @rdname analyze
setMethod("analyze", signature("data.frame"),
          function(data, estimator, data_distribution, use_full_twoarm_sampling_distribution, design, sigma, exact){
            if (missing(data))
              stop("data argument may not be ommited.")
            if (is(data_distribution, "Student") && !missing(sigma)){
              warning("data_distribution was set to Normal because sigma was specified.")
              data_distribution <- Normal(two_armed = data_distribution@two_armed)
            }
            sdata <- summarize_data(data)
            if (missing(sigma)){
              sigma <- sqrt(sdata$svar)
              if (is(data_distribution, "Normal")){
                warning("data_distribution is of class Normal but sigma was not specified. Estimating sigma from the data.")
              }
            }
            arglist <- c(sdata, design = design, sigma = sigma, two_armed = data_distribution@two_armed)
            if (!is.list(estimator)){
              estimator <- list(estimator)
            }
            results <- lapply(estimator, .analzye,
                              data_distribution = data_distribution,
                              use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                              design = design,
                              sigma = sigma,
                              exact = exact,
                              arglist = arglist)
            return(
              Results(data = data,
                      summary_data = sdata,
                      data_distribution = data_distribution,
                      design = design,
                      sigma = sigma,
                      results = results)
              )
  })
.analzye <- function(estimator, data_distribution, use_full_twoarm_sampling_distribution, design, sigma, exact, arglist){
  stagewise_estimators <- get_stagewise_estimators(estimator = estimator,
                                                   data_distribution =  data_distribution,
                                                   use_full_twoarm_sampling_distribution = use_full_twoarm_sampling_distribution,
                                                   design = design, sigma = sigma, exact = exact)
  ret <- list(estimator = estimator)
  if (arglist$n_stages==2L){
    if (is(estimator, "IntervalEstimator")) {
      resl <- do.call(stagewise_estimators[[3L]], arglist)
      resr <- do.call(stagewise_estimators[[4L]], arglist)
      res <- list(lower = resl, upper = resr)
    } else {
      res <- list(do.call(stagewise_estimators[[2L]], arglist))
    }
    ret[["stage2"]] <- res
  } else{
    if (is(estimator, "IntervalEstimator")) {
      resl <- do.call(stagewise_estimators[[1L]], arglist)
      resr <- do.call(stagewise_estimators[[2L]], arglist)
      res <- list(lower = resl, upper = resr)
    } else {
      res <- list(do.call(stagewise_estimators[[1L]], arglist))
    }
    ret[["stage1"]] <- res
  }
  ret
}
summarize_data <- function(data, reference_value = 0, endpoint_name = "endpoint", group_name = "group", stage_name = "stage"){
  if (!is.data.frame(data)){
    stop("data needs to be a data.frame.")
  }
  for (nm in c(endpoint_name, group_name, stage_name)){
    if (! (nm %in% names(data)) ){
      stop(sprintf("no variable with the name %s is present in data.", nm))
    }
  }
  names(data)[names(data) %in% c(endpoint_name, group_name, stage_name)] <- c(endpoint_name, group_name, stage_name)
  if (!is.factor(data$group))
    stop("group needs to be a factor.")
  if (length(levels(data$group)) > 2L )
    stop("group needs to have 1 or 2 levels.")
  if (!is.numeric(data$endpoint))
    stop("endpoints needs to be numeric")
  if (!all(data$stage %in% c(1L,2L)))
    stop("stage may only contain the numbers 1 and 2.")
  n_stages <- length(unique(data$stage))
  n_groups <- length(levels(data$group))
  data$group <- as.numeric(data$group)
  data_s1 <- data[data$stage==1L,]
  data_s1_g1 <- data_s1[data_s1$group==1L,]
  n_s1_g1 <- sum(!is.na(data_s1_g1$endpoint))
  n1 <- n_s1_g1
  ret <- data.frame(n_stages = n_stages,
                    n_groups = n_groups,
                    n_s1_g1 = n_s1_g1)
  smean1T <- mean(data_s1_g1$endpoint, na.rm=TRUE)
  if (n_groups == 2L){
    data_s1_g2 <- data_s1[data_s1$group==2L,]
    n_s1_g2 <- sum(!is.na(data_s1_g2$endpoint))
    n1 <- n1 + n_s1_g2
    smean1C <- mean(data_s1_g2$endpoint, na.rm=TRUE)
    smean1 <- smean1T - smean1C
    ret <- cbind(ret, n_s1_g2 = n_s1_g2, n1 = n1, smean1T = smean1T, smean1C = smean1C, smean1 = smean1)
  } else{
    smean1 = smean1T - reference_value
    ret <- cbind(reference_value = reference_value, ret, n1 = n1, smean1T = smean1T, smean1 = smean1)
  }
  if (n_groups == 2L){
    svar1 <- ((n_s1_g1 - 1L) * var(data_s1_g1$endpoint, na.rm=TRUE) + (n_s1_g2 - 1L) * var(data_s1_g2$endpoint, na.rm=TRUE)) / (n_s1_g1 + n_s1_g2 - 2L)
  } else{
    svar1 <- var(data_s1_g1$endpoint, na.rm=TRUE)
  }
  ret <- cbind(ret, svar1 = svar1)
  if (n_stages==2L){
    data_s2 <- data[data$stage==2L,]
    data_s2_g1 <- data_s2[data_s2$group==1L,]
    n_s2_g1 <- sum(!is.na(data_s2_g1$endpoint))
    n2 <-  n_s2_g1
    ret <- cbind(ret, n_s2_g1 = n_s2_g1)
    smean2T <- mean(data_s2_g1$endpoint, na.rm=TRUE)
    if (n_groups == 2L){
      data_s2_g2 <- data_s1[data_s2$group==2L,]
      n_s2_g2 <- sum(!is.na(data_s2_g2$endpoint))
      n2 <- n2 + n_s2_g2
      smean2C <- mean(data_s2_g2$endpoint, na.rm=TRUE)
      smean2 <- smean2T - smean2C
      ret <- cbind(ret, n_s2_g2 = n_s2_g2, n2 = n2, smean2T = smean2T, smean2C = smean2C, smean2 = smean2)
    } else{
      smean2 = smean2T - reference_value
      ret <- cbind(ret, n_s2_g1 = n_s2_g1, n2 = n2, smean2T = smean2T, smean2 = smean2)
    }
    if (n_groups == 2L){
      svar2 <- ((n_s2_g1 - 1L) * var(data_s2_g1$endpoint, na.rm=TRUE) + (n_s2_g2 - 1L) * var(data_s2_g2$endpoint, na.rm=TRUE)) / (n_s2_g1 + n_s2_g2 - 2L)
    } else{
      svar2 <- var(data_s2_g1$endpoint, na.rm=TRUE)
    }
    ret <- cbind(ret, svar2 = svar2)
  }
  if (n_groups == 2L){
    data_g1 <- data[data$group==1L,]
    n_g1 <- sum(!is.na(data_g1$endpoint))
    data_g2 <- data[data$group==2L,]
    n_g2 <- sum(!is.na(data_g2$endpoint))
    svar <- (var(data_g1$endpoint) * (n_g1 - 1L) + var(data_g2$endpoint) * (n_g2 - 1L) ) / (n_g1 + n_g2 - 2L)
    data_s1_g1 <- data_s1[data_s1$group==1L,]
  } else{
    data_g1 <- data[data$group==1L,]
    n_g1 <- sum(!is.na(data_g1$endpoint))
    svar <- var(data_g1$endpoint)
  }
  ret <- cbind(ret, svar = svar)
  ret
}
