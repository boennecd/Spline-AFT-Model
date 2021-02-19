# source the file with the estimation code
source("fit-saft.R")

# prepare colon data set
colon_use <- within(
  colon[complete.cases(colon) & colon$etype==2, ], {
    time <- time/365.25
    age_c <- scale(age)
  })

mf <- model.frame(time ~ rx + sex + age_c + obstruct + perfor + adhere +
                    factor(differ) + factor(extent) + surg + node4,
                  colon_use)

X <- model.matrix(terms(mf), mf)

y <- model.response(mf)
event <- colon_use$status

# get Gauss–Legendre quadrature node and weights
gl_dat <- get_gl_rule(50)

# estimate model using B-splines
system.time(
  res_bs <- saft_fit(y = y, X = X, event = event, n_knots = 2, gl_dat = gl_dat,
                     basis_type = "bs"))
res_bs$coefs
res_bs$logLik

# we can also use C++ which is faster
system.time(
  res_bs_cpp <- saft_fit_cpp(y = y, X = X, event = event, n_knots = 2,
                             gl_dat = gl_dat, basis_type = "bs"))

all.equal(res_bs[c("coefs", "par")], res_bs_cpp[c("coefs", "par")])

# estimate the model using natural cubic splines
system.time(
  res_ns <- saft_fit(y = y, X = X, event = event, n_knots = 4, gl_dat = gl_dat,
                     basis_type = "ns"))
res_ns$coefs
res_ns$logLik

# again, we can also use C++
system.time(
  res_ns_cpp <- saft_fit_cpp(y = y, X = X, event = event, n_knots = 4,
                             gl_dat = gl_dat, basis_type = "ns"))

all.equal(res_ns[c("coefs", "par")], res_ns_cpp[c("coefs", "par")])

# plot the two estimated splines
jpeg("saft_fit-spline.jpg", width = 600, height = 400)
ys <- seq(0, 2 * max(y), length.out = 1000)
par(mar = c(5, 5, 1, 1))
with(res_bs, plot(
  ys, eval_basis(ys) %*% gamma, type = "l", ylim = c(-6, -2.75),
  main = "", xlab = "Time", ylab = "Spline", bty = "l", yaxs = "i", xaxs = "i"))
with(res_ns, lines(ys, eval_basis(ys) %*% gamma, lty = 2))
grid()
dev.off()

# plot the estimated hazards
jpeg("saft_fit-hazard.jpg", width = 600, height = 400)
ys <- seq(0, 2 * max(y), length.out = 1000)
par(mar = c(5, 5, 1, 1))
with(res_bs, plot(
  ys, exp(eval_basis(ys) %*% gamma), type = "l", ylim = c(0, .08),
  main = "", xlab = "Time", ylab = "Baseline hazard", bty = "l", yaxs = "i",
  xaxs = "i"))
with(res_ns, lines(ys, exp(eval_basis(ys) %*% gamma), lty = 2))
grid()
dev.off()

# compare with Weibull model. Should be the same as Pang et al. (2021).
survreg_res <- survreg(Surv(y, event) ~ X - 1)
coef(survreg_res)
logLik(survreg_res)

# check that we get the same with integrate
system.time(
  res_ns_int <- saft_fit(y = y, X = X, event = event, n_knots = 4,
                         gl_dat = gl_dat, use_integrate = TRUE,
                         basis_type = "ns", check_grads = FALSE))
rbind(res_ns$coefs, res_ns_int$coefs)
c(`Gauss–Legendre quadrature` = res_ns$logLik,
  `R's integrate` = res_ns_int$logLik)
