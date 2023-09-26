pad_middle <- function(left, right, maxlen = 80L){
  len <- maxlen - nchar(left) - nchar(right)
  if (len > 0L){
    return(paste0(left, paste0(rep(" ", len), collapse=""), right))
  } else{
    return(paste0(left, " ", right))
  }
}

#' @importFrom utils capture.output
setMethod("toString", signature("EstimatorScoreResult"),
          function(x, maxlen = 80L, ...){
            lines <- list()
            left <- "Design:"
            right <- toString(x@design)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            left <- "Data Distribution:"
            right <- toString(x@data_distribution)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            left <- "Estimator:"
            right <- x@estimator@label
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            left <- "Assumed sigma:"
            right <- x@sigma
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            left <- "Assumed mu:"
            right <- paste(format(x@mu, ...), collapse = " ")
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            lines[[length(lines)+1L]] <- paste0("Results:\n")
            for (i in seq_along(x@results)){
              res <- x@results[[i]]
              nm <- names(x@results)[[i]]
              left <- paste0(" ", nm, collapse="")
              right <- paste(format(res, ...), collapse = " ")
              lines[[length(lines)+1L]] <- paste0(pad_middle(paste0(left, ":"), right, maxlen = maxlen), "\n")
            }
            return(unlist(lines))
          })
setMethod("show", signature("EstimatorScoreResult"), \(object) cat(c(toString(object), "\n"), sep=""))

setMethod("toString", signature("Statistic"),
          function(x, ...){
            return(paste0(x@label))
          })
setMethod("show", signature("Statistic"), \(object)cat(c(toString(object), "\n"), sep=""))

#' @importFrom utils capture.output
setMethod("toString", signature("Results"),
          function(x, ...) {
            lines <- list()
            left <- "Design:"
            right <- toString(x@design)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Data Distribution:"
            right <- toString(x@data_distribution)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Observed number of stages:"
            right <- x@summary_data$n_stages
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            if (x@data_distribution@two_armed){
              left <- "Observed n1 (group 1)"
              nval1 <- x@summary_data$n_s1_g1
              right <- format(nval1)
              lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

              left <- "Observed n1 (group 2)"
              nval1 <- x@summary_data$n_s1_g2
              right <- format(nval1)
              lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")
            }
            left <- "Observed n1 (total)"
            nval1 <- x@summary_data$n1
            right <- format(nval1)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            test_str <- if (is(x@data_distribution, "Normal")) "Z1" else "T1"
            if (x@summary_data$n_groups == 2L)
              n1 <- c(x@summary_data$n_s1_g1, x@summary_data$n_s1_g2)
            else
              n1 <- x@summary_data$n1

            test_val <-
              if (is(x@data_distribution, "Normal"))
                z_test(x@summary_data$smean1,
                       n1,
                       x@sigma,
                       x@data_distribution@two_armed)
            else
              t_test(x@summary_data$smean1,
                     x@summary_data$svar1,
                     n1,
                     x@data_distribution@two_armed)
            left <- test_str
            right <- format(test_val, digits=3)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Interim decision:"
            if (test_val > x@design@c1e) {
              right <- "reject null (early efficacy stop)"
            } else if (test_val< x@design@c1f) {
              right <- "accept null (early futility stop)"
            } else {
              right <- "continue to second stage"
            }
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Calculated n2(Z1) (per group)"
            nval2 <- .n2_extrapol(x@design, test_val)
            right <- format(nval2)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Calculated c2(Z1)"
            cval2 <- .c2_extrapol(x@design, test_val)
            right <- format(cval2, digits=3)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            if (x@summary_data$n_stages==2L){
              if (x@data_distribution@two_armed){
                left <- "Observed n2 (group 1)"
                nval2 <- x@summary_data$n_s2_g1
                right <- format(nval2)
                lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

                left <- "Observed n2 (group 2)"
                nval2 <- x@summary_data$n_s2_g2
                right <- format(nval2)
                lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")
              }
              left <- "Observed n2 (in total)"
              nval2 <- x@summary_data$n2
              right <- format(nval2)
              lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

              test_str <- if (is(x@data_distribution, "Normal")) "Z2" else "T2"
              if (x@summary_data$n_groups == 2L)
                n2 <- c(x@summary_data$n_s2_g1, x@summary_data$n_s2_g2)
              else
                n2 <- x@summary_data$n2
              test_val2 <-
                if (is(x@data_distribution, "Normal"))
                  z_test(x@summary_data$smean2,
                         n2,
                         x@sigma,
                         x@data_distribution@two_armed)
              else
                t_test(x@summary_data$smean2,
                       x@summary_data$svar2,
                       n2,
                       x@data_distribution@two_armed)
              left <- test_str
              right <- format(test_val2, digits=3)
              lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

              left <- "Final test decision:"
              if (test_val2 > cval2) {
                right <- "reject null"
              } else {
                right <- "accept null"
              }
              lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")
            }
            print_header <- TRUE
            for (res in x@results){
              if ("stage1" %in% names(res)){
                if (print_header){
                  lines[[length(lines)+1L]] <- paste0("\n")
                  lines[[length(lines)+1L]] <- paste0("Stage 1 results:\n")
                  print_header <- FALSE
                }
                left <- paste0(" ", res$statistic@label, collapse="")
                right <- format(res$stage1, ...)
                lines[[length(lines)+1L]] <- paste0(pad_middle(paste0(left, ":"), right), "\n")
              }
            }
            print_header <- TRUE
            if (x@summary_data$n_stages==2L){
              for (res in x@results){
                if ("stage2" %in% names(res)){
                  if (print_header){
                    lines[[length(lines)+1L]] <- paste0("\n")
                    lines[[length(lines)+1L]] <- paste0("Stage 2 results:\n")
                    print_header <- FALSE
                  }
                  left <- paste0(" ", res$statistic@label, collapse="")
                  if (length(res$stage2) > 1L){
                    right <- paste0("[",paste0(format(res$stage2, ...), collapse = ", "), "]")
                  } else {
                    right <- format(res$stage2, ...)
                  }
                  lines[[length(lines)+1L]] <- paste0(pad_middle(paste0(left, ":"), right), "\n")
                }
              }
            }
            return(unlist(lines))
          })
setMethod("show", signature("Results"), \(object) cat(c(toString(object), "\n"), sep = ""))
setMethod("toString", signature("DataDistribution"),
          function(x, ...) {
            sprintf("%s<%s>", class(x), if(x@two_armed) "two-armed" else "single-armed")
})
setMethod("toString", signature("EstimatorScore"),
          function(x, ...) {
            x@label
          })
setMethod("toString", signature("TwoStageDesign"),
          function(x, ...) {
            if (!is.null(attr(x, "label")))
              return(attr(x, "label"))
            # This is partially taken from the adoptr implementation of design2str
            n2_piv <- seq(x@c1f, x@c1e, length.out = length(x@n2_pivots))
            n2range <- range(.n2_extrapol(x, n2_piv))
            sprintf(
              "%s<n1=%i;%.1f<=x1<=%.1f:n2=%s>",
              class(x)[1], n1(x, round=TRUE),
              x@c1f, x@c1e,
              if (diff(round(n2range)) == 0) sprintf("%i", round(n2range)[1]) else paste(round(n2range), collapse = '-')
            )
          })
setMethod("toString", signature("TwoStageDesignWithCache"),
          function(x, ...) {
            toString(forget_cache(x), ...)
          })

### Remove once adoptr is back on CRAN ###
setMethod("show", signature("TwoStageDesign"),
          function(object) {
            cat(toString(object))
          })
### end remove ###

setGeneric("toTeX", \(x, ...) standardGeneric("toTeX"))

#' @importFrom latex2exp TeX
setMethod("toTeX", signature("ANY"),
          function(x, ...) {
          toString(x, ...)
          })
#' @importFrom latex2exp TeX
setMethod("toTeX", signature("NeymanPearsonOrderingPValue"),
          function(x, ...) {
            str <- sprintf("Neyman-Pearson test ordering ($\\mu_0 = %.1f, \\mu_1 = %.1f$)", x@mu0, x@mu1)
            str
          })


format.EstimatorScoreResultList <- function(x, ...) rep("<EstimatorScoreResult>", length(x))
`[.EstimatorScoreResultList` <- function(x, i){
  class(x) <- class(x)[class(x)!="EstimatorScoreResultList"]
  x <- x[i]
  class(x) <- c("EstimatorScoreResultList", class(x))
  x
}

setMethod("toString", "NormalPrior",
          \(x, ...) {
            sprintf("NormalPrior<mu=%s;sigma=%s>", format(x@mu), format(x@sigma))
          })
setMethod("show", signature("NormalPrior"), \(object)cat(c(toString(object), "\n"), sep=""))

setMethod("toString", "UniformPrior",
          \(x, ...) {
            sprintf("UniformPrior<min=%s;max=%s>", format(x@min), format(x@max))
          })
setMethod("show", signature("UniformPrior"), \(object)cat(c(toString(object), "\n"), sep=""))

