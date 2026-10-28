# Runs CmdStan's gradient test for `model` at one point and returns the
# model and finite-difference gradients (on the unconstrained scale), plus
# flags for the two ways initialisation can fail.
#
# Uses `$diagnose()` so cmdstanr handles paths and the runtime environment
# on every platform. CmdStan exits with an error when any finite-difference
# error exceeds `error`, so it is set high here and the comparison is left
# to the tests. A failed initialisation still errors, and is reported
# through the flags, read from the captured CmdStan output.
stan_gradient_at <- function(model, data, init) {
  fit <- NULL
  out <- utils::capture.output(
    fit <- tryCatch( # nolint: implicit_assignment_linter.
      model$diagnose(data = data, init = list(init), error = 1e10),
      error = function(e) NULL
    )
  )
  grads <- if (is.null(fit)) NULL else fit$gradients()
  list(
    rejected = any(grepl("Rejecting initial value", out, fixed = TRUE)),
    gradient_not_finite = any(grepl(
      "Gradient evaluated at the initial value", out,
      fixed = TRUE
    )),
    gradient = if (is.null(grads)) numeric(0) else grads$model,
    finite_diff = if (is.null(grads)) numeric(0) else grads$finite_diff
  )
}

# Returns the log density and its gradient on the unconstrained scale at
# `init`. `model` must be compiled with `compile_model_methods = TRUE`.
stan_log_prob_and_grad <- function(model, data, init) {
  fit <- suppressMessages(suppressWarnings(model$sample(
    data = data, init = list(init), chains = 1, iter_sampling = 1,
    fixed_param = TRUE, refresh = 0, show_messages = FALSE
  )))
  fit$init_model_methods(verbose = FALSE)
  upars <- fit$unconstrain_variables(init)
  list(log_prob = fit$log_prob(upars), grad = fit$grad_log_prob(upars))
}
