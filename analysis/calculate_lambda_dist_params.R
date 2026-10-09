# Table 1. Distributions of the model transition rates
#
# Reads each parameter set's input file and returns, for every rate with a
# probabilistic uncertainty model, the lognormal parameters used by
# sample_lambda() (Supplement A.2.2):
#
#   Median survival time x with standard error se (years), point_stat = "median":
#     X ~ LogN(mu, sigma^2), mu = log(x), sigma = se / x
#     lambda = log(2) / X ~ LogN(log(log(2)) - mu, sigma^2)
#   Incidence rate x with SD of the log rate se (1/year), point_stat = "ir":
#     lambda ~ LogN(log(x), se^2)
#
# Input:  Data/input_params_<set>.csv
# Output: Fig/Table1_lambda_moments.csv

devtools::load_all()
library(data.table)

param_files <- c(
  "NETTER&RADIANT" = "Data/input_params_NETTER&RADIANT.csv",
  "COMPETE" = "Data/input_params_COMPETE.csv"
)

# Table 1 symbols, in table order
param_symbols <- c(
  lambda_everolimus_12 = "lambda_12^E",
  lambda_prrt_12 = "lambda_12^P",
  lambda_cabozantinib_23 = "lambda_23^C",
  lambda_46 = "lambda_46",
  lambda_dot5 = "lambda_.5",
  lambda_57 = "lambda_57"
)

#' Lognormal parameters of lambda, as sampled by sample_lambda()
lambda_lognormal_moments <- function(est, se, point_stat) {
  stopifnot(point_stat %in% c("median", "ir"))
  meanlog <- fifelse(point_stat == "median", log(log(2)) - log(est), log(est))
  sdlog <- fifelse(point_stat == "median", se / est, se)
  list(meanlog = meanlog, varlog = sdlog^2, sdlog = sdlog)
}

moments <- rbindlist(lapply(names(param_files), function(set_name) {
  param_input <- as.data.table(read.csv(param_files[[set_name]]))
  param_input <- param_input[param %in% names(param_symbols)]
  param_input[, c("meanlog", "varlog", "sdlog") := lambda_lognormal_moments(est, se, point_stat)]
  param_input[, `:=`(
    parameter_set = set_name,
    symbol = param_symbols[param],
    # Summaries of lambda itself (1/year)
    lambda_median = exp(meanlog),
    lambda_mean = exp(meanlog + varlog / 2),
    lambda_lower95 = qlnorm(0.025, meanlog, sdlog),
    lambda_upper95 = qlnorm(0.975, meanlog, sdlog),
    table1 = sprintf("LogN(%.3f, %.3f)", meanlog, varlog)
  )]
  param_input[, .(parameter_set, param, symbol, point_stat, est, se, meanlog, varlog, sdlog,
                  lambda_median, lambda_mean, lambda_lower95, lambda_upper95, table1, source)]
}))
moments[, param := factor(param, levels = names(param_symbols))]
setorder(moments, param, parameter_set)

write.csv(moments, file.path("Fig", "Table1_lambda_dist_params.csv"), row.names = FALSE)

# Table 1 entries, LogN(mu, sigma^2), one column per parameter set
print(dcast(moments, param + symbol ~ parameter_set, value.var = "table1"))

# Check against the model: log of the rates drawn by sample_parameters()
set.seed(20260926)
check <- rbindlist(lapply(names(param_files), function(set_name) {
  param_input <- as.data.table(read.csv(param_files[[set_name]]))
  draws <- sample_parameters(N = 1e6, param_input = param_input)
  data.table(
    parameter_set = set_name,
    param = names(param_symbols),
    meanlog_sampled = vapply(names(param_symbols), function(p) mean(log(draws[[p]])), numeric(1)),
    varlog_sampled = vapply(names(param_symbols), function(p) var(log(draws[[p]])), numeric(1))
  )
}))
check <- merge(moments[, .(parameter_set, param = as.character(param), meanlog, varlog)], check,
               by = c("parameter_set", "param"))
message(sprintf(
  "Check against sample_parameters(): max |meanlog diff| = %.4f, max |varlog diff| = %.4f",
  max(abs(check$meanlog - check$meanlog_sampled)),
  max(abs(check$varlog - check$varlog_sampled))
))
