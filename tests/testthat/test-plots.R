plot_chains <- function() {
  set.seed(11)
  X <- lapply(1:2, function(v) {
    m <- matrix(rnorm(40 * 2, rep(c(0, 3), each = 20)), 40, 2)
    rownames(m) <- 1:40
    m
  })
  withr_quiet <- options(mdir.quiet = TRUE)
  on.exit(options(withr_quiet))
  fitMDI(X, 3, R = 200, thin = 5, types = c("G", "MVN"), K = c(4, 4), verbose = FALSE)
}

built <- function(p) {
  expect_s3_class(p, "ggplot")
  ggplot2::ggplot_build(p)
}

test_that("plot() of fits gives a ggplot for every type", {
  ch <- plot_chains()
  for (type in c("trace", "density", "rhat", "psm", "fusion")) {
    expect_error(built(plot(ch, type = type)), NA, info = type)
  }
  expect_error(built(plot(ch[[1]])), NA)
  expect_error(built(plot(ch, type = "psm", view = 2)), NA)
  expect_error(built(plot(ch, type = "trace", pars = c("phi[1,2]", "mass[1]"))), NA)
})

test_that("trace plots omit the initial state and mark the burn in", {
  ch <- plot_chains()
  p <- plot(ch, type = "trace")
  expect_true(all(p$data$iteration > 0))
  expect_equal(sort(unique(p$data$iteration)), seq(5, 200, by = 5))
  expect_equal(nlevels(p$data$chain), 3)
  vl <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  expect_equal(vl[[1]]$data$xintercept, 100)
  expect_equal(levels(p$data$quantity),
               c("complete_likelihood", "joint_likelihood", "mass[1]", "occupied components[1]",
                 "mass[2]", "occupied components[2]", "phi[1,2]"))
})

test_that("density plots use only the retained draws", {
  ch <- plot_chains()
  p <- plot(ch, type = "density", burn = 100)
  # 200 / 5 + 1 = 41 saved; the initial state and the first 20 are dropped
  expect_equal(max(table(p$data$chain, p$data$quantity)), 20)
  expect_true(all(p$data$iteration > 100))
})

test_that("processed chains are plotted at their true iterations", {
  ch <- plot_chains()
  pr <- processMCMCChains(ch, burn = 100)
  p <- plot(pr, type = "trace")
  expect_equal(min(p$data$iteration), 105)
  expect_equal(max(p$data$iteration), 200)
  expect_error(built(plot(pr, type = "psm")), NA)
  vl <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  expect_length(vl, 0)
})

test_that("psm and fusion plots hold the expected numbers", {
  ch <- plot_chains()
  p <- plot(ch, type = "psm", view = 1)
  expect_equal(nrow(p$data), 40^2)
  expect_true(all(p$data$probability >= 0 & p$data$probability <= 1))
  expect_equal(sum(p$data$probability[p$data$item_x == p$data$item_y]), 40)
  f <- plot(ch, type = "fusion")
  expect_equal(nrow(f$data), 40)
  # the sorted values are those of calcFusionProbability on the pooled retained draws
  idx <- floor(100 / 5) + 2
  pooled <- sort(colMeans(do.call(rbind, lapply(unclass(ch), function(c) c$allocations[idx:41, , 1])) ==
    do.call(rbind, lapply(unclass(ch), function(c) c$allocations[idx:41, , 2]))))
  expect_equal(f$data$probability, pooled)
})

test_that("plot() fails clearly on bad requests", {
  ch <- plot_chains()
  expect_error(plot(ch, pars = "nope"), "Unknown quantity")
  expect_error(plot(ch, type = "psm", view = 3), "view")
  one <- callMDI(list(matrix(rnorm(40), 20, 2)), R = 40, thin = 5, types = "G", K = 3, check_prior = FALSE)
  expect_s3_class(plot(one), "ggplot")
  expect_error(plot(one, type = "fusion"), "two views")
})

test_that("convergence and smc objects have plot methods", {
  ch <- plot_chains()
  conv <- assessConvergence(ch, burn = 100)
  p <- plot(conv)
  expect_error(built(p), NA)
  expect_equal(sum(p$data$series == "Rhat"), nrow(conv))
  set.seed(2)
  X <- list(matrix(rnorm(60, rep(c(0, 3), each = 15)), 30, 2))
  rownames(X[[1]]) <- 1:30
  s <- smcMDI(X, "G", K = 3, n_particles = 30)
  expect_error(built(plot(s)), NA)
  expect_equal(nrow(plot(s)$data), 3 * nrow(s$trace))
})
