# Helpers for testing the C++ versions of Stan functions (see
# pcd_stan_header()) against their Stan bodies.

# Compiles `model_code` (a model without a functions block) twice: with the
# Stan bodies, and with the C++ functions declared and the header included.
# Returns a named list of the two CmdStanModel objects.
compile_cpp_and_stan <- function(model_code, name) {
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
  )
  build <- function(cpp) {
    functions <- pcd_load_stan_functions(wrap_in_block = TRUE, cpp = cpp)
    model_dir <- file.path(tempdir(), paste0(name, if (cpp) "cpp" else "stan"))
    dir.create(model_dir, showWarnings = FALSE)
    path <- file.path(model_dir, paste0(name, ".stan"))
    code <- paste0(functions, "\n", model_code)
    # Rewriting an unchanged file would make cmdstanr recompile the model
    code_lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
    if (!file.exists(path) || !identical(readLines(path), code_lines)) {
      writeLines(code, path)
    }
    model_args <- list(path, dir = model_dir)
    if (cpp) {
      model_args$user_header <- pcd_stan_header()
    }
    suppressMessages(suppressWarnings(
      do.call(cmdstanr::cmdstan_model, model_args)
    ))
  }
  list(stan = build(FALSE), cpp = build(TRUE))
}

# Skips unless CmdStan's log_prob method can be run directly. It needs
# CmdStan 2.31 or later, and on Windows the binary needs the TBB library on
# the path, which cmdstanr only sets up for its own runs.
skip_if_no_log_prob <- function() {
  testthat::skip_on_os("windows")
  testthat::skip_if(
    cmdstanr::cmdstan_version(error_on_NA = FALSE) < "2.31"
  )
}

# Log density and its gradient at parameters `upars` (a named list), from
# CmdStan's log_prob method at full precision. The probe models have only
# unconstrained parameters, so these are on the unconstrained scale.
stan_log_prob <- function(model, data, upars) {
  run_dir <- tempfile("logprob")
  dir.create(run_dir)
  data_file <- file.path(run_dir, "data.json")
  pars_file <- file.path(run_dir, "pars.json")
  out_file <- file.path(run_dir, "out.csv")
  cmdstanr::write_stan_json(data, data_file)
  cmdstanr::write_stan_json(upars, pars_file)
  run_log <- suppressWarnings(system2(
    model$exe_file(),
    c(
      "log_prob", paste0("constrained_params=", pars_file),
      "data", paste0("file=", data_file),
      "output", paste0("file=", out_file), "sig_figs=18"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(run_log, "status"))) {
    stop(paste(run_log, collapse = "\n"), call. = FALSE)
  }
  out <- utils::read.csv(out_file, comment.char = "#")
  list(lp = out$lp__[1], grad = unlist(out[1, -1], use.names = FALSE))
}

# Generated quantities named `name` at parameters `init`, at full precision.
stan_gq <- function(model, data, init, name) {
  fit <- model$sample(
    data = data, init = list(init), fixed_param = TRUE, chains = 1,
    iter_sampling = 1, sig_figs = 18, refresh = 0, show_messages = FALSE,
    show_exceptions = FALSE
  )
  as.vector(fit$draws(name, format = "draws_matrix"))
}
