# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# installs branches to benchmark
touchstone::branch_install()

touchstone::benchmark_run(
  expr_before_benchmark = { source("touchstone/setup.R") },
  <benchmark_expr> = { <expr_to_benchmark> },
  n = 5
)

# create artifacts used downstream in the GitHub Action.
touchstone::benchmark_analyze()
