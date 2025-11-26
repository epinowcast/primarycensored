if (requireNamespace("spelling", quietly = TRUE)) {
  # Skip spelling checks if _R_CHECK_SPELLING_ is set to false
  # (used in CI for non-release R versions)
  run_spelling <- Sys.getenv("_R_CHECK_SPELLING_", "true") != "false"
  if (run_spelling) {
    spelling::spell_check_test(
      vignettes = TRUE,
      error = TRUE,
      skip_on_cran = TRUE,
      lang = "en-GB"
    )
  }
}
