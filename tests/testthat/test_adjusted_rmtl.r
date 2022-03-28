
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
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5789, 0.4279))
})

test_that("rmtl surv, no boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, interpolation="linear")
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5957, 0.4401))
})

test_that("rmtl surv, no boot but use_boot=TRUE", {
  adj_rmtl <- suppressWarnings(adjusted_rmtl(adj_no_boot, to=1.1,
                                             use_boot=TRUE))
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5789, 0.4279))
})

test_that("rmtl surv, with boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5789, 0.4279))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.0511, 0.0587))
  expect_equal(as.vector(adj_rmtl$n_boot), c(7, 5))
})

test_that("rmtl surv, with boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, use_boot=TRUE, interpolation="linear")
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5957, 0.4401))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.0501, 0.0586))
  expect_equal(as.vector(adj_rmtl$n_boot), c(7, 5))
})

test_that("rmtl surv, no boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5435, 0.4191))
})

test_that("rmtl surv, with boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.5435, 0.4191))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.0420, 0.0525))
  expect_equal(as.vector(adj_rmtl$n_boot), c(7, 5))
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
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1856, 0.1599))
})

test_that("rmtl cif, no boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, interpolation="linear")
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1918, 0.1618))
})

test_that("rmtl cif, with boot", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1856, 0.1599))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.1185, 0.0775))
  expect_equal(as.vector(adj_rmtl$n_boot), c(10, 10))
})

test_that("rmtl cif, with boot, linear", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, use_boot=TRUE,
                            interpolation="linear")
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1918, 0.1618))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.1204, 0.0784))
  expect_equal(as.vector(adj_rmtl$n_boot), c(10, 10))
})

test_that("rmtl cif, no boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1650, 0.1536))
})

test_that("rmtl cif, with boot, using from", {
  adj_rmtl <- adjusted_rmtl(adj, to=1.1, from=0.3, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmtl$auc, 4)), c(0.1650, 0.1536))
  expect_equal(as.vector(round(adj_rmtl$auc_se, 4)), c(0.1104, 0.0664))
  expect_equal(as.vector(adj_rmtl$n_boot), c(10, 10))
})
