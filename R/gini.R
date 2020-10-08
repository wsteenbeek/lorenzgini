#' @import boot
#' @import DescTools
#' @import ggplot2
NULL

#' Calculates the (generalized) Gini coefficient
#'
#' This function calculates the Gini coefficient. It builds on the
#' `gini` function from package `DescTools`, and extends this
#' function to the generalized Gini.
#'
#' See 'Bernasco, W. and W. Steenbeek (2017). More places than crimes:
#' Implications for evaluating the law of crime concentration at place.
#' Journal of Quantitative Criminology. https://doi.org/10.1007/s10940-016-9324-7'
#' for more information.
#'
#' @param mydata     a vector of N units, with counts of events per unit. The
#' vector should contain at least non-negative elements. The result will be NA,
#' if the vector contains negative elements.
#' @param generalized  logical. Should the generalized gini be calculated?
#' (default = TRUE)
#' @param unbiased     logical. Should the unbiased Gini be calculated?
#' (default = TRUE)
#' @param conf.level   confidence level for the confidence interval, restricted
#' to lie between 0 and 1. If set, bootstrap confidence intervals are
#' calculated. If set to NA (default) no confidence intervals are returned.
#' @param R            number of bootstrap replicates. Usually this will be a
#' single positive integer. For importance resampling, some resamples may use
#' one set of weights and others use a different set of weights. In this case
#' R would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights. This is ignored if no confidence
#' intervals are to be calculated.
#' @param type         character string representing the type of interval
#' required. The value should be one out of the c("norm","basic", "stud",
#' "perc" or "bca"). This argument is ignored if no confidence intervals are
#' to be calculated.
#' @param plot         logical. Should a kernel density plot of the distribution
#' of R Ginis be printed? (default = FALSE). This argument is ignored if no
#' confidence intervals are to be calculated.
#' @return             If conf.level is set to NA then the result will be a
#' single numeric value and if a conf.level is provided, a named numeric vector
#' with 3 elements: gini (Gini coefficient), lwr.ci	(lower bound of the
#' confidence interval), upr.ci	(upper bound of the confidence interval)
#' @examples
#'
#' crimes <- c(10,8,2,0,0,rep(0,5))
#' crimes
#' length(crimes) # 10 units of analysis
#' sum(crimes) # 20 crimes in total
#' # Standard Gini:
#' gini(crimes, generalized = FALSE, unbiased = FALSE) # 0.78
#' # Is the same the generalized Gini, as Ncrimes >= Nunits:
#' gini(crimes, generalized = TRUE, unbiased = FALSE) # 0.78
#'
#' crimes <- c(5,2,1,0,0,rep(0,5))
#' crimes
#' length(crimes) # 10 units of analysis
#' sum(crimes) # 8 crimes in total
#' # Standard Gini:
#' gini(crimes, generalized = FALSE, unbiased = FALSE) # 0.8
#' # Is not the same the generalized Gini, as Ncrimes < Nunits:
#' gini(crimes, generalized = TRUE, unbiased = FALSE) # 0.75
#'
#' # Example of unbiased standard Gini:
#' crimes <- c(10, rep(0,9))
#' crimes
#' length(crimes) # 10 units of analysis
#' sum(crimes) # 10 crimes in total
#' gini(crimes, generalized = FALSE, unbiased = FALSE) # 0.9
#' # All crime is concentrated in one unit, but the Gini does not equal 1.
#'
#' crimes <- c(10, rep(0,99))
#' crimes
#' gini(crimes, generalized = FALSE, unbiased = FALSE) # 0.99
#'
#' # As the number of units with 0 crimes inceases, the Gini asymptotically
#' # approaches 1. This is adjusted in the 'unbiased' Gini, the default in `DescTools`:
#' crimes <- c(10, rep(0,9))
#' gini(crimes, generalized = FALSE, unbiased = TRUE) # 1
#' crimes <- c(10, rep(0,99))
#' gini(crimes, generalized = FALSE, unbiased = TRUE) # 1
#'
#' # generate crime events from zero-inflated poisson disribution
#' set.seed(8722)
#' crimes <- ifelse(rbinom(1000, size = 1, prob = .1) > 0, 0, rpois(1000, lambda = 1))
#'
#' # frequency of crimes
#' table(crimes)
#'
#' # Bootstrap-based confidence intervals:
#' set.seed(347)
#' gini(crimes, conf.level = .95)
#' # Also output kernel density plot of generated Ginis, with estimated Gini
#' # indicated as solid red line, and the confidence intervals delineated
#' # by dashed red lines:
#' set.seed(347)
#' gini(crimes, conf.level = .95, plot = TRUE)
#' # increasing number of simulations
#' set.seed(347)
#' gini(crimes, conf.level = .95, R = 20000, plot = TRUE)
#'
#' @export
gini <- function(mydata, generalized = TRUE, unbiased = TRUE, conf.level = NA, R = 1000, type = "bca", plot = FALSE){

  # from DescTools::Gini
  # cast to numeric, as else sum(mydata * 1:n) might overflow for integers
  # http://stackoverflow.com/questions/39579029/integer-overflow-error-using-gini-function-of-package-desctools
  mydata <- as.numeric(mydata)

  # from DescTools::Gini
  if (any(is.na(mydata)) || any(mydata < 0)) return(NA_real_)

  mygini <- function(x, generalized = TRUE, unbiased = TRUE){

    Nunits <- length(x)
    Ncrimes <- sum(x)

    # if we don't want the Gini for sparse data, simply
    # run the DescTools::Gini function
    if(!generalized){
      output <- DescTools::Gini(x, unbiased = unbiased)
    }

    # for Gini for sparse data:
    if(generalized){

      # return biased generalized Gini
      if(!unbiased){
        Gini.biased <- DescTools::Gini(x, unbiased = unbiased)
        output <- max(1, Nunits/Ncrimes) * (Gini.biased - 1) + 1
      }

      # return unbiased Gini for sparse data
      if(unbiased){
        Gini.biased <- DescTools::Gini(x, unbiased = FALSE)
        Gini.gen.biased <- max(1, Nunits/Ncrimes) * (Gini.biased - 1) + 1
        Gini.gen.unbiased <- Gini.gen.biased * min(Nunits, Ncrimes) / (min(Nunits, Ncrimes) - 1)
        output <- Gini.gen.unbiased
      }
    }

    return(output)

  }

  if(is.na(conf.level)){
    res <- mygini(mydata, generalized = generalized, unbiased = unbiased)
  } else {
    # adjusted bootstrap percentile (BCa) interval
    boot.gini <- boot::boot(mydata, function(x, d) mygini(x[d], generalized = generalized, unbiased = unbiased), R = R)
    ci <- boot::boot.ci(boot.gini, conf = conf.level, type = type)

    if(plot) {
      # how much should x-axis be shown left from mean and right from mean?
      dev <- max(abs(range(boot.gini$t) - boot.gini$t0))
      xlimits <- c(boot.gini$t0 - dev, boot.gini$t0 + dev)

      myplot <- ggplot2::ggplot(data = data.frame(id = seq_along(boot.gini$t), gini = boot.gini$t), aes(x = gini)) +
        geom_density(fill = "#4271AE", colour = "#1F3552",
                     alpha = 0.6) +
        scale_x_continuous(name = paste0("N = ", length(boot.gini$t))) +
        scale_y_continuous(name = "Density") +
        geom_vline(xintercept = boot.gini$t0, linetype="solid", color = "red") +
        geom_vline(xintercept = ci[[4]][4], linetype="dashed", color = "red") +
        geom_vline(xintercept = ci[[4]][5], linetype="dashed", color = "red") +
        expand_limits(x = xlimits) +
        ggtitle("Kernel density plot of distribution of Ginis")

      print(myplot)
    }

    res <- c(gini = boot.gini$t0, lwr.ci = ci[[4]][4], upr.ci = ci[[4]][5])
  }

  return(res)
}
