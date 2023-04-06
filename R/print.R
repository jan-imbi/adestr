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
            right <- substr(.tmp <- capture.output(print(x@design)), 1L, nchar(.tmp)-1L)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right, maxlen = maxlen), "\n")

            left <- "Data Distribution:"
            right <- substr(.tmp <- capture.output(print(x@data_distribution)), 1L, nchar(.tmp)-1L)
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

setMethod("toString", signature("Estimator"),
          function(x, ...){
            return(paste0(x@label))
          })
setMethod("show", signature("Estimator"), \(object)cat(c(toString(object), "\n"), sep=""))

#' @importFrom utils capture.output
setMethod("toString", signature("Results"),
          function(x, ...) {
            lines <- list()
            left <- "Design:"
            right <- substr(.tmp <- capture.output(print(x@design)), 1L, nchar(.tmp)-1L)
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            left <- "Data Distribution:"
            right <- substr(.tmp <- capture.output(print(x@data_distribution)), 1L, nchar(.tmp)-1L)
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

            left <- "Actual number of stages:"
            right <- x@summary_data$n_stages
            lines[[length(lines)+1L]] <- paste0(pad_middle(left, right), "\n")

            print_header <- TRUE
            for (res in x@results){
              if ("stage1" %in% names(res)){
                if (print_header){
                  lines[[length(lines)+1L]] <- paste0("\n")
                  lines[[length(lines)+1L]] <- paste0("Stage 1 results:\n")
                  print_header <- FALSE
                }
                left <- paste0(" ", res$estimator@label, collapse="")
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
                  left <- paste0(" ", res$estimator@label, collapse="")
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

#' @importFrom utils capture.output
setMethod("toString", signature("DataDistribution"),
          function(x, ...) {
            str <- capture.output(print(x))
            substr(str, 1, nchar(str)-1L)
})
setMethod("toString", signature("EstimatorScore"),
          function(x, ...) {
            x@label
          })
#' @importFrom utils capture.output
setMethod("toString", signature("TwoStageDesign"),
          function(x, ...) {
            if (!is.null(attr(x, "label")))
              return(attr(x, "label"))
            str <- capture.output(print(x))
            substr(str, 1, nchar(str)-1L)
          })

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

#' @export
format.EstimatorScoreResultList <- function(x, ...) rep("<EstimatorScoreResult>", length(x))

#' @export
`[.EstimatorScoreResultList` <- function(x, i){
  class(x) <- class(x)[class(x)!="EstimatorScoreResultList"]
  x <- x[i]
  class(x) <- c("EstimatorScoreResultList", class(x))
  x
}



