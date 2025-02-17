# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Simulate leapfrog model
#'
#' @param demp list of demographic input parameters (TODO: document)
#' @param projp list of HIV projection parameters (TODO: document)
#' @param hiv_strat stratification of HIV population, either "full"
#'   (default; single-year ages) or "coarse" (aggregated age groups). 
#' @param hiv_steps_per_year number of Euler integration steps per year
#'   for HIV progression; default 10.
#'
#' @details
#' The first year of `sx`, `asfr`, `srb`, and `netmig` is not used. This is assumed
#' to apply to the base year population (consistent with Spectrum).
#'
#' @export
#' 
leapfrogR <- function(demp, projp, hiv_strat = "full", hiv_steps_per_year = 10L) {
    .Call(`_leapfrog_leapfrogR`, demp, projp, hiv_strat, hiv_steps_per_year)
}

