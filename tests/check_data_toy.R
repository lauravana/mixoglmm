tol <- 1e-4
# devtools::install_github("lauravana/mixoglmm")
library(mixoglmm)
library(optimx)
data("data_toy", package = "mixoglmm")

######### Test cor_general ##########
## cor_general
formula <- (Be2 + Be1 + y1 + y2 + y3 ~ 1 + X1 + X2)
system.time(
fit <- mixoglmm(formula,
                families = list(
                  Be1 = binomial(link="probit"),
                  Be2 = binomial(link="probit"),
                  y1 = gaussian(),
                  y2 = gaussian(),
                  y3 = gaussian()),
                data = data_toy,
                constraints.beta = list("(Intercept)" =
                                     cbind(c(1,1,1,1,1),
                                           c(0,1,0,0,0),
                                           c(0,0,1,0,0),
                                           c(0,0,0,1,0),
                                           c(0,0,0,0,1)),
                                   X1 = cbind(c(1,1,1,1,1)),
                                   X2 = cbind(c(1,1,1,1,1))),
                cor_struct_gauss = cor_general(~ 1),
                control = mixoglmm.control(solver = "nlminb"),
                na.action = "na.pass")
)

summary(fit)

p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

options(digits = 22)
ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(-0.0148211422205021, 0.0658555379858056, -0.84682781285835, 0.168396865029936, 1.06806874085768, -0.953346745451055, 0.445402687630554, 0.931367272712356, 0.956898907763614, 1.05984932467885, 1.06234349259962, 1.13755789952594, 0.948380641180452, 0.54522461347693, 0.411477527885525, 0.373503876518674, 0.924098693958423, 2.01334905278068),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.081448617788922, 0.0988458750995877, 0.0719298898156584, 0.0812656742698099, 0.114254366005471, 0.0518345026426737, 0.047956637718494, 0.139771492748076, 0.141764242006868, 0.0913037856440705, 0.107737972865219, 0.171403917432732, 0.184297520945764, 0.162779637223033, 0.0799824338027637, 0.236349520217856, 0.111384019354783, 0.101896420374603),
                           tol = tol, check.attributes = FALSE))

## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:10], c(-0.468462989439877, -0.182413011204749, -2.02827297987693, -2.20107649932213, 0.284808999755551, 0.0835657596957073, 1.51969881481497, -1.24086011355168, -1.55930573018537, -0.124845711512893),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.468761979746234, -0.182452170713698, -2.02836483520674, -2.18752870444045, 0.285526019496651, 0.0827293593327543, 1.51973502355592, -1.24118899651708, -1.5596704185102, -0.126555479460619),
                        tol = tol, check.attributes = FALSE))

######### Test cor_equi ##########
formula <- (Be1 + y1 + y2 + y3 ~ 1 + X1 + X2)
system.time(
  fit <- mixoglmm(formula,
                  families = list(
                    Be1 = binomial(link="probit"),
                    y1 = gaussian(),
                    y2 = gaussian(),
                    y3 = gaussian()),
                  data = data_toy,
                  constraints.beta = list("(Intercept)" =
                                            cbind(c(1,1,1,1),
                                                  c(0,1,0,0),
                                                  c(0,0,1,0),
                                                  c(0,0,0,1)),
                                          X1 = cbind(c(1,1,1,1)),
                                          X2 = cbind(c(1,1,1,1))),
                  constraints.lambda = list(cbind(c(1,1,1,1))),
                  control = mixoglmm.control(solver = "nlminb"),
                  cor_struct_gauss = cor_equi())
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$gauss.corr, ss$gauss.stddev)
mixoglmm:::check(all.equal(params[, 1], c(0.0518289069276356, -0.91034189834756, 0.104882391488756, 1.0045498838863, -0.959879197007644, 0.467626408106982, 1.10474655755588, 0.30902792841799, 0.165219063487385, 0.620966169773287, 1.89688070153668),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.0866136895502943, 0.0711627086371993, 0.07603114090254, 0.11048008501755, 0.0542946480948415, 0.0503941401459024, 0.0378149671036537, 0.090261697655027, 0.0732986786925112, 0.0328661380492024, 0.0726101277995174),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:10], c(-0.342249300764135, -0.198829717174722, -1.52713525946909, -2.50265409159475, -0.207361843519006, 0.382829696651032, 1.55369549993144, -1.01316475399982, -0.653746967539385, -0.2185919748955),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.342326525202454, -0.198747928569311, -1.52718463275267, -2.50263792035243, -0.207443200401388, 0.382913651977901, 1.55370601950216, -1.01323384371224, -0.653830965119969, -0.218656838797524),
                           tol = tol, check.attributes = FALSE))

######### Test cor_ident ##########
system.time(
  fit <- mixoglmm(formula,
                  families = list(
                    Be1 = binomial(link="probit"),
                    y1 = gaussian(),
                    y2 = gaussian(),
                    y3 = gaussian()),
                  data = data_toy,
                  constraints.beta = list("(Intercept)" =
                                            cbind(c(1,1,1,1),
                                                  c(0,1,0,0),
                                                  c(0,0,1,0),
                                                  c(0,0,0,1)),
                                          X1 = cbind(c(1,1,1,1)),
                                          X2 = cbind(c(1,1,1,1))),
                  control = mixoglmm.control(solver = "nlminb"),
                  constraints.lambda = list(cbind(c(1,1,1,1))),
                  cor_struct_gauss = cor_ident())
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)
mixoglmm:::check(all.equal(params[, 1], c(0.0509029665688094, -0.909083435525103, 0.106140867726725, 1.00580906593832, -0.96184682595998, 0.469390184162426, 1.12347211416172, 2.69352076297138e-05, 0.582060995456782, 1.83823555114918),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.0869573475088768, 0.070602966909314, 0.0752487145315025, 0.10836511830494, 0.0546248010080437, 0.0507080874981021, 0.035527308495373, 0.0303076549929852, 0.0184063848064434, 0.058130107871574),
                           tol = tol, check.attributes = FALSE))

## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:10], c(-0.352290934571053, -0.224051282553091, -1.55158315420583, -2.56699220460811, -0.203974632003134, 0.430353045305748, 1.5629581697742, -1.01519838511296, -0.635858126572852, -0.224648605670073),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.352290934571054, -0.224051282553091, -1.55158315420583, -2.56699220460811, -0.203974632003134, 0.430353045305748, 1.5629581697742, -1.01519838511296, -0.635858126572852, -0.224648605670074),
                           tol = tol, check.attributes = FALSE))


######### Other link functions for binomial ########

  ###### Logit #####

formula2 <- (Be1 + y1 ~ X1)
fit_logit <- mixoglmm(formula2,
                families = list(
                  Be1 = binomial(link="logit"),
                  y1 = gaussian()),
                data = data_toy,
                constraints.beta = list("(Intercept)" = cbind(c(1,1)),
                                   X1 = cbind(c(1,1))),
                constraints.lambda = list(cbind(c(1,1))),
                control = mixoglmm.control(solver = "nlminb"),
                cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_logit, method = "conditional means")

p_mode <- extract_ranef(fit_logit, method = "conditional modes")

ss <- summary(fit_logit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

# paste0(params[, 1], collapse = ", ")
mixoglmm:::check(all.equal(params[, 1], c(-0.918514266214344, -0.925170501787154, 1.10269398324232, 9.18094415337041e-05),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.0544851665208486, 0.0589647043553358, 0.017435127050066, 0.12346383351973),
                           tol = tol, check.attributes = FALSE))
# paste0(p_mode[1:5], collapse = ", ")
mixoglmm:::check(all.equal(p_mode[1:5], c(-1.13291075522549, 0.325366027782728, -2.00297573617848, -3.12360963053506, -1.02490309995612),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-1.13291075522549, 0.325366027782728, -2.00297573617848, -3.12360963053506, -1.02490309995612),
                           tol = tol, check.attributes = FALSE))

###### Cloglog #####
fit_cloglog <- mixoglmm(formula2,
                      families = list(
                        Be1 = binomial(link="cloglog"),
                        y1 = gaussian()),
                      data = data_toy,
                      constraints.beta = list("(Intercept)" = cbind(c(1,1)),
                                              X1 = cbind(c(1,1))),
                      constraints.lambda = list(cbind(c(1,1))),
                      control = mixoglmm.control(solver = "nlminb"),
                      cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_cloglog, method = "conditional means")

p_mode <- extract_ranef(fit_cloglog, method = "conditional modes")

ss <- summary(fit_cloglog)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(-0.905563201906023, -0.930438375291139, 1.06157654308154, 0.299834880150429),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.0555410293227458, 0.0590948822730252, 0.0378369818709493, 0.120561240232432),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mode[1:5], c(-1.10785480778146, 0.332931764769181, -1.9429705642293, -2.94673507587798, -1.00573767886636),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-1.10818122319653, 0.332608912270392, -1.94315566472562, -2.94683835767591, -1.00612419682856),
                           tol = tol, check.attributes = FALSE))



###### Cauchit #####
fit_cauchit <- mixoglmm((Be1 + y3 ~ X1 + X2),
                        families = list(
                          Be1 = binomial(link="cauchit"),
                          y3 = gaussian()),
                        data = data_toy,
                        constraints.beta = list("(Intercept)" = cbind(c(1,1),
                                                                      c(0,1)),
                                                X1 = cbind(c(1,1)),
                                                X2 = cbind(c(1,1))),
                        constraints.lambda = list(cbind(c(1,1))),
                        control = mixoglmm.control(solver = "nlminb"),
                        cor_struct_gauss = cor_general(~ 1))

p_mean <- extract_ranef(fit_cauchit, method = "conditional means")

p_mode <- extract_ranef(fit_cauchit, method = "conditional modes")

ss <- summary(fit_cauchit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.0213560469166707, 1.05776039372046, -1.15430283708733, 0.560250392079567, 0.730279839373291, 1.33161313195184),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.114090673257652, 0.136949784826862, 0.098988603781712, 0.086699839575671, 0.0868329726349795, 0.0514517209558047),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-0.37781367090812, -0.383972959593408, -0.943813209113344, -0.700001461244408, -0.187361077881396),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-0.394862632372445, -0.364158330243734, -0.968389686127352, -0.707326920607082, -0.207449727634569),
                           tol = tol, check.attributes = FALSE))



# #######################
# ## Factor covariates ##
# #######################
# data(data_toy)
#
#
#
# data_toy$X4 <- factor(sample(c("d", "e", "f"), n, replace = TRUE))
#
# formula <-  Be1 + y3 ~ 1 + X1 + X2 + X3 + X4
#
# fit <- mixoglmm(formula,
#                 families = list(
#                   Be1 = binomial(link="probit"),
#                   y3 = gaussian()),
#                 data = data_toy,#,
#                 contrasts =  list(X3 = function(x)
#                   contr.treatment(nlevels(data_toy$X3), base = 2),
#                   X4 = "contr.sum"),
#                 # constraints = list("(Intercept)" =
#                 #                      cbind(c(1,1),
#                 #                            c(0,1)),
#                 #                    X1 = (cbind(c(1,1))),
#                 #                    X2 = (cbind(c(1,1))),
#                 #                    X31 = (cbind(c(1,1))),
#                 #                    X33 = (cbind(c(1,1))),
#                 #                    X41= (cbind(c(1,1))),
#                 #                    X42 =(cbind(c(1,1)))),
#
#                 cor_struct_gauss = cor_general(~ 1))
# fit$gradient
# fit$hessian
# fit$se
#
# head(fit$x)
# str(fit$mf)
# ################################################################
# fit1 <- mixoglmm(Be1 + y1 ~ 1 + X1 + X2,
#                 #     na.action = "na.pass",
#                 families = list(
#                   Be1 = binomial(link="probit"),
#                   # Be2 = binomial(link="probit"),
#                  # Po1 = poisson(),
#                   y1 = gaussian()),
#                 #  y2 = gaussian(),
#                 #  y3 = gaussian()),
#                 data = data_toy,
#                 cor_struct_gauss = cor_general(~ 1))
#
# ####
# library(mcglm)
#
# Z0_ex4 <- mc_id(data_toy)
# fit.joint.log <- mcglm(linear_pred = c(Be1 ~ X1 + X2, y3 ~ X1 + X2),
#                        matrix_pred = list(Z0_ex4,Z0_ex4),
#                        link = c("probit", "identity"),
#                        data = data_toy)
#
# #########################
# library(mcglm)
# data("soya", package = "mcglm")
# form.grain <- grain ~ block + water * pot
# form.seed <- seeds ~ block + water * pot
# soya$viablepeasP <- soya$viablepeas / soya$totalpeas
# form.peas <- viablepeasP ~ block + water * pot
# Z0_ex4 <- mc_id(soya)
# fit.grain <- mcglm(linear_pred = c(form.grain), matrix_pred = list(Z0_ex4), data = soya)
# fit.seed <- mcglm(linear_pred = c(form.seed), matrix_pred = list(Z0_ex4),
#                   link = c("log"), variance = c("tweedie"), power_fixed = TRUE, data = soya)
# fit.peas <- mcglm(linear_pred = c(form.peas), matrix_pred = list(Z0_ex4),
#                   link = "logit", variance = "binomialP", Ntrial = list(soya$totalpeas),
#                   data = soya)
#
# fit.joint <- mcglm(linear_pred = c(form.grain, form.seed, form.peas),
#
#                     matrix_pred = list(Z0_ex4, Z0_ex4, Z0_ex4), link = c("identity",
#                    "log", "logit"), variance = c("constant", "tweedie", "binomialP"),
#                    Ntrial = list(NULL, NULL, soya$totalpeas), data = soya)
# summary(fit.joint, verbose = TRUE, print = "Correlation")
#
# mm <- model.matrix( ~ block + water * pot, data = soya)
# colnames(mm)
# #constr <- rep(list(cbind(c(1,1,1))), NCOL(mm))
# #names(constr) <- colnames(mm)
# fit.soya <- mixoglmm(grain + seeds + viablepeas ~ block + water * pot,#+ offset(totalpeas),
#                      #     na.action = "na.pass",
#                      families = list(
#                        grain  = gaussian(),
#                        # Be2 = binomial(link="probit"),
#                        seeds = poisson(),
#                        viablepeas = binomial(link="logit")),
#                      Ntrials = list(NULL, NULL, soya$totalpeas),
#
#                     # offset = list(0.2*soya$totalpeas, NULL, NULL),
#                      #  y2 = gaussian(),
#                      #  y3 = gaussian()),
#                      data = soya,
#                     # constraints = constr,
#                      cor_struct_gauss = cor_general(~ 1))
#
# extract_ranef(fit.soya, "conditional modes")
#
# gof(fit.peas)
# gof(fit.grain)
# gof(fit.seed)

######### Missing values ########

set.seed(12345)
data_toy_NA <- data_toy
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "y1"] <- NA
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "y2"] <- NA
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "Be1"] <- NA

## cor_general
formula <- (Be1 + y1 + y2 + y3 ~ 1 + X1 + X2)
system.time(
  fit_NA <- mixoglmm(formula, na.action = na.pass,
                  families = list(
                    Be1 = binomial(link="probit"),
                    y1 = gaussian(),
                    y2 = gaussian(),
                    y3 = gaussian()),
                  data = data_toy_NA,
                  constraints.beta = list("(Intercept)" =
                                            cbind(c(1,1,1,1),
                                                  c(0,1,0,0),
                                                  c(0,0,1,0),
                                                  c(0,0,0,1)),
                                          X1 = cbind(c(1,1,1,1)),
                                          X2 = cbind(c(1,1,1,1))),
                  constraints.lambda = list(cbind(c(1,1,1,1))),
                  cor_struct_gauss = cor_general(~ 1),
                  control = mixoglmm.control(solver = "nlminb"))
)

p_mean_NA <- extract_ranef(fit_NA, method = "conditional means")

p_mode_NA <- extract_ranef(fit_NA, method = "conditional modes")

options(digits = 22)
ss <- summary(fit_NA)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.0154682071008168, -0.901533872195918, 0.115816837002234, 1.0218370206747, -0.919572126743135, 0.446605235923827, 0.932973483039596, 0.836947936660481, 0.526977741355371, 0.463455153562708, 0.6003456425333, 0.981122200097407, 2.08157406796045),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.0832528767851219, 0.0798171184590386, 0.0879919701643304, 0.120213453371418, 0.0516550527554206, 0.0477194216387664, 0.0590807633643879, 0.0289024094722624, 0.056908517549467, 0.0501545671621834, 0.0811766782970644, 0.0710316398457183, 0.07874072564092),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode_NA[1:5], c(-0.383646713249743, -0.0551746308961111, -1.59624494499529, -1.66635271170662, -0.300998141256397),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean_NA[1:5], c(-0.389663714341792, -0.049042132970619, -1.59624494499529, -1.66395870860669, -0.307139801512674),
                           tol = tol, check.attributes = FALSE))

######### Missing values for all normal responses ########
set.seed(12345)
data_toy_NA <- data_toy
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "y1"] <- NA
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "y2"] <- NA
data_toy_NA[sample(1:nrow(data_toy_NA), 50), "Be1"] <- NA
data_toy_NA[1:2, c("y1", "y2", "y3")] <- NA

## cor_general
formula <- (Be1 + y1 + y2 + y3 ~ 1 + X1 + X2)
system.time(
  fit_NA <- mixoglmm(formula, na.action = na.pass,
                     families = list(
                       Be1 = binomial(link="probit"),
                       y1 = gaussian(),
                       y2 = gaussian(),
                       y3 = gaussian()),
                     data = data_toy_NA,
                     constraints.beta = list("(Intercept)" =
                                               cbind(c(1,1,1,1),
                                                     c(0,1,0,0),
                                                     c(0,0,1,0),
                                                     c(0,0,0,1)),
                                             X1 = cbind(c(1,1,1,1)),
                                             X2 = cbind(c(1,1,1,1))),
                     constraints.lambda = list(cbind(c(1,1,1,1))),
                     cor_struct_gauss = cor_general(~ 1),
                     control = mixoglmm.control(solver = "nlminb"))
)

p_mean_NA <- extract_ranef(fit_NA, method = "conditional means")

p_mode_NA <- extract_ranef(fit_NA, method = "conditional modes")

options(digits = 22)
ss <- summary(fit_NA)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.0154699510660333, -0.900471632077768, 0.116198263259389, 1.02626618845146, -0.920191281908805, 0.446748486377923, 0.934838734324552, 0.836760430237387, 0.526352642417378, 0.463357733186428, 0.601184345511543, 0.98286647053191, 2.08440266855998),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.0832704675295902, 0.0798662713348333, 0.0880997380731996, 0.120469458420921, 0.0518114359465066, 0.0479747738559247, 0.0592574867743312, 0.0289921656042577, 0.0570821135827522, 0.0502563174669352, 0.0813986705477454, 0.0712512852151604, 0.0790053897556807),
                           tol = tol, check.attributes = FALSE))

## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode_NA[1:5], c(-0.0573582928052127, 0.0365629660613572, -1.59403516681913, -1.66380339367687, -0.301061097848),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean_NA[1:5], c(-0.0636471758921174, 0.0426048325714321, -1.59403516681913, -1.66139948875571, -0.307229905919011),
                           tol = tol, check.attributes = FALSE))
