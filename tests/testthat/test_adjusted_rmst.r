
suppressMessages(requireNamespace("survival"))

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

test_that("rmst 2 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adj, to=1.1)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5211, 0.6721))
})

test_that("rmst 2 treatments, no boot, linear", {
  adj_rmst <- adjusted_rmst(adj, to=1.1, interpolation="linear")
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.50430, 0.6599))
})

test_that("rmst 2 treatments, no boot but use_boot=TRUE", {
  adj_rmst <- suppressWarnings(adjusted_rmst(adj_no_boot, to=1.1,
                                             use_boot=TRUE))
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5211, 0.6721))
})

test_that("rmst 2 treatments, with boot", {
  adj_rmst <- adjusted_rmst(adj, to=1.1, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5211, 0.6721))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0511, 0.0587))
  expect_equal(as.vector(adj_rmst$n_boot), c(7, 5))
})

test_that("rmst 2 treatments, with boot, linear", {
  adj_rmst <- adjusted_rmst(adj, to=1.1, use_boot=TRUE,
                            interpolation="linear")
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5043, 0.6599))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0501, 0.0586))
  expect_equal(as.vector(adj_rmst$n_boot), c(7, 5))
})

test_that("rmst 2 treatments, no boot, using from", {
  adj_rmst <- adjusted_rmst(adj, to=1.1, from=0.3)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.2565, 0.3809))
})

test_that("rmst 2 treatments, with boot, using from", {
  adj_rmst <- adjusted_rmst(adj, to=1.1, from=0.3, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.2565, 0.3809))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0420, 0.0525))
  expect_equal(as.vector(adj_rmst$n_boot), c(7, 5))
})

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=sum(sim_dat$group==1),
                                          replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

test_that("rmst 3 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adj, to=1)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5115, 0.7021, 0.6293))
})

test_that("rmst 3 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adj, to=1)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5115, 0.7021, 0.6293))
})

test_that("rmst 3 treatments, with boot, linear", {
  adj_rmst <- adjusted_rmst(adj, to=1, use_boot=TRUE, interpolation="linear")
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.4958, 0.6867, 0.6221))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0691, 0.0628, 0.0703))
  expect_equal(as.vector(adj_rmst$n_boot), c(10, 6, 8))
})

test_that("rmst 3 treatments, no boot, using from", {
  adj_rmst <- adjusted_rmst(adj, to=1, from=0.3)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.2470, 0.4021, 0.3441))
})

test_that("rmst 3 treatments, with boot, using from", {
  adj_rmst <- adjusted_rmst(adj, to=1, from=0.3, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.2470, 0.4021, 0.3441))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0582, 0.0639, 0.0599))
  expect_equal(as.vector(adj_rmst$n_boot), c(10, 6, 8))
})
