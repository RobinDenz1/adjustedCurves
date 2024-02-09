
suppressMessages(requireNamespace("survival"))

### using single-event data

set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

adj_no_boot <- adj
adj_no_boot$boot_data <- NULL
adj_no_boot$boot_adjsurv <- NULL

test_that("rmtl surv, no boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5789, 0.4279))
})

test_that("rmtl surv, no boot, diff", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, difference=TRUE)
  expect_equal(round(adj_rmtl$diff, 4), 0.151)
})

test_that("rmtl surv, no boot, ratio", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, ratio=TRUE)
  expect_equal(round(adj_rmtl$ratio, 4), 1.3529)
})

test_that("rmtl surv, no boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, interpolation="linear")
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5957, 0.4401))
})

test_that("rmtl surv, no boot but conf_int=TRUE", {
  adj_rmtl <- suppressWarnings(adjusted_rmtl(adj_no_boot, to=1.1,
                                             conf_int=TRUE))
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5789, 0.4279))
})

test_that("rmtl surv, with boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5789, 0.4279))
  expect_equal(round(adj_rmtl$se, 3), c(0.078, 0.068))
  expect_equal(adj_rmtl$n_boot, c(10, 7))
})

test_that("rmtl surv, with boot, diff", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE, difference=TRUE)
  expect_equal(round(adj_rmtl$diff, 4), 0.151)
  expect_equal(round(adj_rmtl$se, 4), 0.1036)
  expect_equal(round(adj_rmtl$ci_lower, 4), -0.0521)
  expect_equal(round(adj_rmtl$ci_upper, 4), 0.3541)
  expect_equal(round(adj_rmtl$p_value, 4), 0.145)
})

test_that("rmtl surv, with boot, ratio", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE, ratio=TRUE)
  expect_equal(round(adj_rmtl$ratio, 4), 1.3529)
  expect_equal(round(adj_rmtl$ci_lower, 4), 0.8984)
  expect_equal(round(adj_rmtl$ci_upper, 4), 2.0969)
  expect_equal(round(adj_rmtl$p_value, 4), 0.145)
})

test_that("rmtl surv, with boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE, interpolation="linear")
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5957, 0.4401))
  expect_equal(round(adj_rmtl$se, 4), c(0.0774, 0.0683))
  expect_equal(adj_rmtl$n_boot, c(10, 7))
})

test_that("rmtl surv, no boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5435, 0.4191))
})

test_that("rmtl surv, with boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3, conf_int=TRUE)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.5435, 0.4191))
  expect_equal(round(adj_rmtl$se, 4), c(0.0663, 0.0580))
  expect_equal(adj_rmtl$n_boot, c(10, 7))
})

### using competing risks data

set.seed(42)
sim_dat <- sim_confounded_crisk(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=10,
                   cause=1)

test_that("rmtl cif, no boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1856, 0.1599))
})

test_that("rmtl cif, no boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, interpolation="linear")
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1918, 0.1618))
})

test_that("rmtl cif, with boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1856, 0.1599))
  expect_equal(round(adj_rmtl$se, 4), c(0.0775, 0.1185))
  expect_equal(adj_rmtl$n_boot, c(10, 10))
})

test_that("rmtl cif, with boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, conf_int=TRUE,
                            interpolation="linear")
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1918, 0.1618))
  expect_equal(round(adj_rmtl$se, 4), c(0.0784, 0.1204))
  expect_equal(adj_rmtl$n_boot, c(10, 10))
})

test_that("rmtl cif, no boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1650, 0.1536))
})

test_that("rmtl cif, with boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3, conf_int=TRUE)
  expect_equal(round(adj_rmtl$rmtl, 4), c(0.1650, 0.1536))
  expect_equal(round(adj_rmtl$se, 4), c(0.0664, 0.1104))
  expect_equal(adj_rmtl$n_boot, c(10, 10))
})
