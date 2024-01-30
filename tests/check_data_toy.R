tol <- 1e-4
####################
# devtools::install_github("lauravana/mixoglmm")
library(mixoglmm)
library(optimx)
data("data_toy", package = "mixoglmm")

## cor_general
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
                cor_struct_gauss = cor_general(~ 1),
                control = mixoglmm.control(solver = "nlminb"),
                na.action = "na.pass")
)

summary(fit)

p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

options(digits = 22)
ss <- summary(fit)
betas <- ss$coefficients
mixoglmm:::check(all.equal(betas[, 1], c(0.158484783477795, -1.00497362867322, -0.0114735251139953, 0.906837185153123, -0.957741880351057, 0.477713149269369),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(betas[, 2], c(0.0708767934490218, 0.0690816934779963, 0.0780015574878195, 0.110764183375751, 0.0278115549330782, 0.0256187432672813),
                           tol = tol, check.attributes = FALSE))
tau <- ss$re.stddev
mixoglmm:::check(all.equal(tau[, 1], c(0.5287064818297940727021),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(tau[, 2], c(0.02598187202672989568053),
                           tol = tol, check.attributes = FALSE))

gamma <- ss$gauss.corr
mixoglmm:::check(all.equal(gamma[, 1], c(0.912465777477577, 0.557376891603283, 0.45042892023527),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(gamma[, 2], c(0.0507491111477424, 0.140651966746755, 0.0732447271450007),
                           tol = tol, check.attributes = FALSE))

omega <- ss$gauss.stddev
mixoglmm:::check(all.equal(omega[,1], c(0.400746553816486, 0.902989385810695, 1.97639376865566),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(omega[,2], c(0.0617308862529334, 0.0694279530110146, 0.0872236299701398),
                           tol = tol, check.attributes = FALSE))
## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:10], c(0.00336105350572743, 0.391727736115974,
                                           -0.721716498058259, -0.766074619316724,
                                           -0.0401553072417684, 0.234808625487268,
                                            0.542031250106634, -0.285321350518882,
                                           -0.416422405299797, -0.36408292291156),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(0.00350347290277018, 0.391965194525711,
                                           -0.722084418004382, -0.76633375648837,
                                           -0.0405232571810766, 0.235173742751026,
                                            0.542280842542217, -0.285678690585022,
                                           -0.416768301204705, -0.364367369820741),
                        tol = tol, check.attributes = FALSE))

## cor_equi
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
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)
mixoglmm:::check(all.equal(params[, 1], c(0.156138650919738, -1.005016745771, -0.0115169813174175, 0.906790157324389, -0.940811562734566, 0.466343692934543, 0.623865010066864, 0.425064831814034, 0.187280559847254, 0.637498362466582, 1.88748737292809),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.0722919534052524, 0.0670845111475402, 0.0724085092926511, 0.107496051343685, 0.0315524688271646, 0.0292504887128358, 0.0155649900950424, 0.0614512681657186, 0.0429320568939804, 0.0327546183449962, 0.0687659552733108),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:10], c(-0.116802808054344, 0.426053650188658, -1.24698845409761, -1.59523917936706, 0.094782598702335, 0.460117032784326, 0.496668696664618, -0.381800665365676, -0.629407140391796, -0.0208800299639348),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.11675925138479, 0.426129008325903, -1.24709452168075, -1.59535559753372, 0.0946649560974065, 0.460235448911945, 0.496753516810754, -0.381917900668997, -0.629524127208963, -0.0209568253587547),
                           tol = tol, check.attributes = FALSE))

## cor_ident
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
mixoglmm:::check(all.equal(params[, 1], c(1.565550788263897208363e-01, -1.005055677165627825431e+00, -1.155608614217086375919e-02,  9.067579782146614197913e-01,
                                          -9.449493971975334583036e-01, 4.674549843605066978824e-01,   6.642033527501332024201e-01,
                                          1.843609848911476311586e-05,  5.629183885736568360159e-01,  1.785009068471818283186e+00),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.07292725470527808284693, 0.06646600503762277456499, 0.07107379297895902170445, 0.10387610209229070079573,
                                          0.03229450514040876213384, 0.02997892096909285722384,  0.01395089666318996210470,
                                          0.01856958157168214079702, 0.01780104259379940156993, 0.05644694413839976970815 ),
                           tol = tol, check.attributes = FALSE))

## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:10], c(-0.26240347611990572085361,  0.50226627730147410666461, -1.37762863206073871147339, -1.78332514418669330069633,
                                           0.06602023442142204712191,  0.50423057625348000065912,  0.50240371047420273598050, -0.50245944550992227561892,
                                           -0.65247293970998165857367,  0.01931721693948460980006  ),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.26240347611990572085361,  0.50226627730147399564231, -1.37762863206073915556260, -1.78332514418669396683015,
                                            0.06602023442142207487748,  0.50423057625348011168143,  0.50240371047420262495820, -0.50245944550992249766352,
                                           -0.65247293970998176959597,  0.01931721693948462367785 ),
                           tol = tol, check.attributes = FALSE))

#######################################
## other link functions for binomial ##
#######################################
## Logit

formula2 <- (Be1 + y1 ~ X1)
fit_logit <- mixoglmm(formula2,
                families = list(
                  Be1 = binomial(link="logit"),
                  y1 = gaussian()),
                data = data_toy,
                constraints.beta = list("(Intercept)" = cbind(c(1,1)),
                                   X1 = cbind(c(1,1))),
                control = mixoglmm.control(solver = "nlminb"),
                cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_logit, method = "conditional means")

p_mode <- extract_ranef(fit_logit, method = "conditional modes")

ss <- summary(fit_logit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

# paste0(params[, 1], collapse = ", ")
mixoglmm:::check(all.equal(params[, 1], c(-0.908585633899121, -0.90842424754779, 0.899860049834889, 0.000270830799918245),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.0362842747538179, 0.0392673999459461, 0.0128036875880352, 0.396948577950672),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-1.11234951443344, 1.11024379553635, -1.84121818499512, -2.34145574116911, -0.83125978690579),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-1.11234951443344, 1.11024379553635, -1.84121818499512, -2.34145574116911, -0.831259786905791),
                           tol = tol, check.attributes = FALSE))

##############
## cloglog
fit_cloglog <- mixoglmm(formula2,
                      families = list(
                        Be1 = binomial(link="cloglog"),
                        y1 = gaussian()),
                      data = data_toy,
                      constraints.beta = list("(Intercept)" = cbind(c(1,1)),
                                              X1 = cbind(c(1,1))),
                      control = mixoglmm.control(solver = "nlminb"),
                      cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_cloglog, method = "conditional means")

p_mode <- extract_ranef(fit_cloglog, method = "conditional modes")

ss <- summary(fit_cloglog)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(-0.8934791679824, -0.916943514097205, 0.840507059325738, 0.322562752501028),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.0381868893225489, 0.0397881886293485, 0.0424152202420524, 0.127722053583628),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-0.896174354654218, 0.966246555527126, -1.62666675155525, -2.08122866117328, -0.751378491961517),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-0.896425961050744, 0.96706952606964, -1.62702932057134, -2.08197922280016, -0.752050050461492),
                           tol = tol, check.attributes = FALSE))



## Cauchit
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

mixoglmm:::check(all.equal(params[, 1], c(0.149226223992418721709, 0.945156020213949732423, -1.241323798981201598224,
                                          0.579684298052147384261,  0.321700803542295687976,
                                          1.388005722277994236080),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.101942986184902661950, 0.133205452009057229157, 0.085435873751396690379,
                                          0.072062241895121767477, 0.051612301832824669656,
                                          0.043391245093964039836),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-0.0815774135194027, 0.0879995177655245, -0.140618826108703, -0.243921134440102, -0.0811467162422348
                                         ),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c( -0.0813012776782469, 0.0883639665123446, -0.141852380642502,
                                           -0.242496184472103, -0.0822578065295178),
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


### Missing values NAs
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

mixoglmm:::check(all.equal(params[, 1], c(0.132776922373239, -0.978883313177678,
                                          0.00588956504176138, 0.931093835352336,
                                          -0.951880150828531, 0.475184594290818,
                                           0.458739249081703, 0.898470356696568,
                                          0.509578939505349, 0.457699468536408,
                                          0.475371743673464, 0.957529445330185, 1.97816425214396),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.0728041125683895, 0.0734415916983707,
                                         0.082949354490307, 0.113231074990345,
                                         0.0287753918023289, 0.0265732094834572,
                                         .0212513284541692, 0.0165462254834365,
                                         0.0473091938811224, 0.0405375893928955, 0.048835468905032, 0.0440211302802085, 0.0661609545313001),
                           tol = tol, check.attributes = FALSE))

tau <- ss$re.stddev
mixoglmm:::check(all.equal(tau[, 1], c(0.4587338756886303881011),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(tau[, 2], c(0.02124770917560831612336 ),
                           tol = tol, check.attributes = FALSE))

gamma <- ss$gauss.corr
mixoglmm:::check(all.equal(gamma[, 1], c(0.8984703566965681620360, 0.5095789395053489378995, 0.4576994685364076387302 ),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(gamma[, 2], c(0.0165462254834365, 0.0473091938811224, 0.0405375893928955 ),
                           tol = tol, check.attributes = FALSE))

omega <- ss$gauss.stddev
mixoglmm:::check(all.equal(omega[,1], c(0.4753717436734636314632, 0.9575294453301854691318, 1.9781642521439577464548 ),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(omega[,2], c(0.04883546890503195209154, 0.04402113028020846613453, 0.06616095453130013681609),
                           tol = tol, check.attributes = FALSE))


## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode_NA[1:5], c(0.0228643722742666,0.325404253151259,-0.524242850627651,-0.586883093910137,0.102140399108836),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean_NA[1:5], c(0.0231618203830685,0.325961459492032,-0.524242850627651,-0.587381861420012,0.101368469936255),
                           tol = tol, check.attributes = FALSE))


### Missing values NAs for all normal responses
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

mixoglmm:::check(all.equal(params[, 1], c(0.132218538111684, -0.978626693128789,
                                          0.0064097339979871, 0.933678624189268,
                                          -0.95045619907437, 0.472146298282196,
                                          0.459974412610806, 0.898826552913782,
                                          0.50696024737938, 0.457434347318497,
                                          0.475033262876167, 0.958934924335357, 1.97503899585749),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.0727961962852914, 0.0734404387232174,
                                         0.0830238494223883, 0.113251093691898,
                                         0.0288472278280454, 0.0267139872440934,
                                         0.0213303869440451, 0.016651922547938, 0.0475522056910673,
                                         0.0406716377851296, 0.0490695083287982, 0.0441898399900465, 0.0662363130748606),
                           tol = tol, check.attributes = FALSE))

## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode_NA[1:5], c(0.108430692069076, 0.00917554058371092, -0.527846230532726, -0.590073477614332, 0.0951891261966428),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean_NA[1:5], c(0.108747283701767, 0.0098635493582526, -0.527846230532727, -0.590572476665559, 0.094417392959938),
                           tol = tol, check.attributes = FALSE))
