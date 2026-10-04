.libPaths(c("/tmp/claude-0/-home-user-MDIr/7d714999-05c1-566d-84c6-cb4710ea967b/scratchpad/lib_new", .libPaths())); suppressMessages(library(mdir)); options(mdir.quiet = TRUE)
source("bimodal_setup.R")
seps <- c(1.0, 1.5, 2.0, 2.5, 3.0)
res <- list()
for (sep in seps) {
  d <- make_data(sep = sep)
  runs <- parallel::mclapply(1:6, function(i) {
    set.seed(300 + i)
    f <- callMDI(list(d$X), R = 8000, thin = 10, types = "MVN", K = 3, check_prior = FALSE, save_parameters = FALSE)
    a <- f$allocations[-(1:100), , 1]
    pat <- apply(a, 1, group_partition, truth = d$truth)
    list(switches = sum(pat[-1] != pat[-length(pat)]), patterns = names(sort(table(pat), decreasing = TRUE))[1:2],
         occ = mean(apply(a, 1, function(l) length(unique(l)))), ll = mean(f$complete_likelihood[-(1:1000)]))
  }, mc.cores = 3)
  cat(sprintf("sep %.1f: plain chains, pattern switches per chain: %s | modal pattern per chain: %s | mean occupied comps %.2f | ll %.0f\n",
    sep, paste(vapply(runs, `[[`, 1, "switches"), collapse = ","), paste(vapply(runs, function(r) r$patterns[1], ""), collapse = " "),
    mean(vapply(runs, `[[`, 1, "occ")), mean(vapply(runs, `[[`, 1, "ll"))))
}
