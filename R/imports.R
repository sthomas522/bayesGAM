#' @importMethodsFrom rstan plot extract show summary 
#' @importFrom cluster clara
#' @importFrom graphics abline
#' @importFrom corrplot corrplot.mixed
# #' @importFrom ggplot2 ggplot aes_string xlab ylab geom_ribbon geom_line facet_wrap theme_bw geom_raster geom_contour labs aes scale_fill_distiller theme_bw element_blank guides
#' @importFrom boot inv.logit
#' @import geometry
#' @import stats
#' @import SemiPar
#' @import ggplot2
NULL

# zzz.R
# .onLoad <- function(libname, pkgname) {
#   if (length(stanmodels) != 0) {
#     modules <- paste0("stan_fit4", names(stanmodels), "_mod")
#     for (m in modules) loadModule(m, what = TRUE)
#   } else {
#     message("No stan programs to compile were found.")
#   }
# }

# TODO: binomial - add weights functionality
# TODO: find binomial example with log link