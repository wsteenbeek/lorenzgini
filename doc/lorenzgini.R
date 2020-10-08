## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = TRUE, eval = TRUE, warning = FALSE, message = FALSE--------------
library("lorenzgini")

## ---- echo = FALSE, fig.show = 'hold', fig.height = 4, fig.width = 4, fig.align = "center", fig.cap="Figure 1: Lorenz curve"----
crimes <- c(22,18,17,12,11,8,6,3,2,1)
lorenz(crimes, prop = FALSE)

## ---- echo = FALSE, fig.show = 'hold', fig.height = 4, fig.width = 4, fig.align = "center", fig.cap="Figure 2: Generalized Lorenz curve, 10 places, 3 + 1 + 1 crimes, G = .78, G’ = .58"----
crimes <- c(3,1,1,rep(0,7))
lorenz(crimes, prop = FALSE)

## ---- eval = FALSE, warning = FALSE, message = FALSE--------------------------
#  library("lorenzgini")

## ---- eval = TRUE, warning = FALSE, message = FALSE, fig.show = 'hold', fig.height = 4, fig.width = 4, fig.align = "center"----
# 10 spatial units of analysis and 100 events
crimes <- c(22,18,17,12,11,8,6,3,2,1)
lorenz(crimes)

## ---- eval = TRUE, warning = FALSE, message = FALSE, fig.show = 'hold', fig.height = 4, fig.width = 4, fig.align = "center"----
# 10 units of analysis and 5 events
crimes <- c(3,1,1,0,0,0,0,0,0,0)
lorenz(crimes)

## ---- eval = TRUE, warning = FALSE, message = FALSE, fig.show = 'hold', fig.height = 4, fig.width = 4, fig.align = "center"----
crimes <- c(3,1,1,rep(0,150))
lorenz(crimes)
lorenz(crimes, rescale = TRUE)

## -----------------------------------------------------------------------------
crimes <- c(10,8,2,0,0,0,0,0,0,0)
length(crimes) # units of analysis
sum(crimes) # crimes

# Standard Gini:
gini(crimes, generalized = FALSE, unbiased = FALSE)

# Is the same the generalized Gini, as Ncrimes >= Nunits:
gini(crimes, generalized = TRUE, unbiased = FALSE)

## -----------------------------------------------------------------------------
crimes <- c(5,2,1,0,0,rep(0,5))
length(crimes) # units of analysis
sum(crimes) # crimes in total

# Standard Gini:
gini(crimes, generalized = FALSE, unbiased = FALSE)

# Is *not* the same the generalized Gini, as Ncrimes < Nunits:
gini(crimes, generalized = TRUE, unbiased = FALSE)

## -----------------------------------------------------------------------------
# Example of unbiased standard Gini:
crimes <- c(10, rep(0,9))
length(crimes) # units of analysis
sum(crimes) # crimes in total

gini(crimes, generalized = FALSE, unbiased = FALSE)
# All crime is concentrated in one unit, but the Gini does not equal 1.

crimes <- c(10, rep(0,99))
crimes
gini(crimes, generalized = FALSE, unbiased = FALSE)

# As the number of units with 0 crimes increases, the Gini asymptotically
# approaches 1. This is adjusted in the 'unbiased' Gini, the default in `DescTools`:
crimes <- c(10, rep(0,9))
gini(crimes, generalized = FALSE, unbiased = TRUE)

crimes <- c(10, rep(0,99))
gini(crimes, generalized = FALSE, unbiased = TRUE)

## ---- fig.height = 4, fig.width = 4, fig.align = "center", message = FALSE, warning = FALSE----
# generate crime events from zero-inflated poisson distribution
set.seed(8722)
crimes <- ifelse(rbinom(1000, size = 1, prob = .1) > 0, 0, rpois(1000, lambda = 1))

# frequency of crimes
table(crimes)

# histogram
ggplot2::qplot(crimes, geom="histogram") 

# Bootstrap-based confidence intervals:
set.seed(347)
gini(crimes, conf.level = .95)

# Also output kernel density plot of generated Ginis, e.g. using 1000 simulations:
gini(crimes, conf.level = .95, R = 1000, plot = TRUE)

# Increasing the number of simulations:
gini(crimes, conf.level = .95, R = 20000, plot = TRUE)

