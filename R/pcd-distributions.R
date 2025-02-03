#' Supported delay distributions
#'
#' @return A data.frame with columns:
#' - name: Distribution name to use with add_name_attribute()
#' - params: Parameter names
#' - aliases: Alternative names/identifiers
#' @export
pcd_distributions <- function() {
  data.frame(
    name = c("lnorm", "gamma", "weibull", "exp"),
    params = c("meanlog,sdlog", "shape,rate", "shape,scale", "rate"),
    aliases = c("lognormal", "gamma", "weibull", "exponential"),
    stan_id = 1:4,
    stringsAsFactors = FALSE
  )
}

#' Supported primary event distributions
#'
#' @return A data.frame with columns:
#' - name: Distribution name to use with add_name_attribute()
#' - params: Parameter names
#' - aliases: Alternative names/identifiers
#' @export
pcd_primary_distributions <- function() {
  data.frame(
    name = c("unif", "expgrowth"),
    params = c("", "r"),
    aliases = c("uniform", "exponential_growth"),
    stan_id = 1:2,
    stringsAsFactors = FALSE
  )
}

#' Get distribution ID by name
#'
#' @param dist_name Distribution name or alias
#' @param type "delay" or "primary"
#'
#' @return Numeric distribution ID
#' @export
pcd_dist_id <- function(dist_name, type = c("delay", "primary")) {
  type <- match.arg(type)
  df <- switch(type,
    delay = pcd_distributions(),
    primary = pcd_primary_distributions()
  )

  match_idx <- which(df$name == dist_name | df$aliases == dist_name)

  if (length(match_idx) == 0) {
    stop("No ", type, " distribution found matching: ", dist_name)
  }

  df$stan_id[match_idx]
}

#' @keywords internal
.suggest_dist_name <- function(input, type = "delay") {
  dist_names <- switch(type,
    delay = c(pcd_distributions()$name, pcd_distributions()$aliases),
    primary = c(
      pcd_primary_distributions()$name,
      pcd_primary_distributions()$aliases
    )
  )

  distances <- utils::adist(input, dist_names)
  min_dist <- min(distances)
  candidates <- dist_names[which(distances == min_dist)]

  if (min_dist <= 2 && length(candidates) > 0) {
    suggestions <- paste0(
      "Did you mean: ",
      toString(unique(candidates)),
      "?"
    )
  } else {
    suggestions <- paste0(
      "Available distributions:",
      toString(unique(dist_names))
    )
  }

  return(suggestions)
}
