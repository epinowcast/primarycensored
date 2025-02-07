#' Supported delay distributions
#'
#' A dataset containing information about the supported delay distributions in
#' primarycensored. Includes both distributions with base R implementations and
#' those only available in Stan. Distributions beyond these are not supported
#' in the stan code but any user functions can be used in the R code.
#'
#' @format A data.frame with 17 rows and 4 columns:
#' \describe{
#'   \item{name}{Distribution name}
#'   \item{pdist}{R distribution function name (e.g. plnorm), NA if there is no
#'   base R implementation}
#'   \item{aliases}{Alternative names/identifiers}
#'   \item{stan_id}{Stan distribution ID used in the stan code}
#' }
#' @family utils
"pcd_distributions"

#' Supported primary event distributions
#'
#' A dataset containing information about the supported primary event
#' distributions in primarycensored. Distributions beyond these are not
#' supported in the stan code but any user functions can be used in the R code.
#'
#' @format A data.frame with 2 rows and 4 columns:
#' \describe{
#'   \item{name}{Distribution name}
#'   \item{dprimary}{R density function name}
#'   \item{aliases}{Alternative names/identifiers}
#'   \item{stan_id}{Stan distribution ID used in the stan code}
#' }
#' @family distribution data
"pcd_primary_distributions"
