# Compare empirical configuration frequencies from several independent chains
# with exact probabilities. Returns z-scores from the between-chain spread.
config_code <- function(mat, K) as.numeric(mat %*% K^(seq_len(ncol(mat)) - 1)) # 0-based labels

freq_by_chain <- function(codes_list, n_codes) {
  vapply(codes_list, function(cd) tabulate(cd + 1, n_codes) / length(cd), numeric(n_codes))
}

compare_to_exact <- function(freq, p_exact, min_p = 2e-3) {
  # freq: n_codes x n_chains
  m <- rowMeans(freq)
  se <- apply(freq, 2, identity) |> apply(1, stats::sd) / sqrt(ncol(freq))
  keep <- p_exact > min_p
  z <- (m[keep] - p_exact[keep]) / se[keep]
  nc <- ncol(freq)
  # between-chain sd estimated from nc chains: the z are t_{nc-1}-like
  tv <- 0.5 * sum(abs(m - p_exact))
  list(tv = tv, n_states = sum(keep), max_abs_z = max(abs(z)),
       mean_z2 = mean(z^2), expected_mean_z2 = (nc - 1) / (nc - 3),
       tv_chain_noise = 0.5 * mean(apply(freq, 2, function(f) sum(abs(f - m)))) / sqrt(nc))
}
