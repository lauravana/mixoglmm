tol <- 1e-6
####################
library(mixoglmm)
data(data_toy)

## cor_general
formula <- (Be1 +  y1 + y2 + y3 ~ 1 + X1 + X2)
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
                cor_struct_gauss = cor_general(~ 1))
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

options(digits = 22)
ss <- summary(fit)
betas <- ss$coefficients
mixoglmm:::check(all.equal(betas[, 1], c(0.15810395787427153346805, -1.00466052811284556334215, -0.01116040877696803171326,
                                          0.90714765855466683586172, -0.95719603771244987644451,  0.47742849573108264715771),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(betas[, 2], c(0.07068215354233951808194, 0.06915542459257395524475, 0.07848438137840953177093,
                                         0.11036586419882268850223, 0.02778705688490125044754, 0.02561253820371892087060),
                           tol = tol, check.attributes = FALSE))
tau <- ss$re.stddev
mixoglmm:::check(all.equal(tau[, 1], c(0.5188263237486250289976),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(tau[, 2], c(0.02179122929320004634079),
                           tol = tol, check.attributes = FALSE))

gamma <- ss$gauss.corr
mixoglmm:::check(all.equal(gamma[, 1], c(0.9321358250970666059487, 0.5261342904515758922557, 0.4405145464523918530375),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(gamma[, 2], c(.02864768955792041521535, 0.05479472138656688723346, 0.04179995102067527584788),
                           tol = tol, check.attributes = FALSE))

omega <- ss$gauss.stddev
mixoglmm:::check(all.equal(omega[,1], c(.4111618152678951365608, 0.9261298856620481245727, 1.9667560746160377416203),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(omega[,2], c(0.05386543183786190319706, 0.04251328757939506347574, 0.06668199753098857929245 ),
                           tol = tol, check.attributes = FALSE))
## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:10], c(-0.03657226644452135633223,  0.41696084576011005484730,
                                           -0.66728045612087527604928, -0.67384220575719233625733,
                                           -0.08249740059866392249965,  0.21076870849528930862427,
                                            0.54600269860613825922968, -0.31662470249903268415537,
                                           -0.37946869284840273328641, -0.40350412012666669880900),
                        tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.03646632621432691762076,  0.41713673889855279464101,
                                           -0.66756018472360567628243, -0.67403006182017211056490,
                                           -0.08277753080010349928664,  0.21104612029620792301721,
                                            0.54619159047551812946608, -0.31689866438335051457287,
                                           -0.37972961720704728749709, -0.40372548667848129344549),
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
                  cor_struct_gauss = cor_equi())
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)
mixoglmm:::check(all.equal(params[, 1], c(0.15613853080842826637209, -1.00501706525141965009595, -0.01151631728663744282104,
                                          0.90679996923632222305400, -0.94081208446587316629461,  0.46634502790876664057862,
                                          1.00000000000000000000000,  0.62386623796274509601290,  0.42506242050030018653217,
                                          0.18727608994427666355698, 0.63749735403715324100204,  1.88748503954429769358114),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.07229196554277539377154, 0.06708447634324424257990, 0.07240848017245393730690,
                                          0.10749596154824477001188, 0.03155247317500815829039, 0.02925049207892149033539,
                                          0.00000000000000000000000, 0.01556510176026645528302, 0.06145285879260477951425,
                                          0.04293345035792595076884, 0.03275504675804622595203, 0.06876636734131037764772 ),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:10], c(-0.11680462800486851593362, 0.42605573208262470519259, -1.24699041281340772258091,
                                           -1.59524213999611497172282, 0.09478470954780673918272,  0.46011961512196714041423,
                                           0.49666922040871896992797, -0.38180231760090155290754,
                                           -0.62940541466724131414168, -0.02087886319661238013201 ),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:10], c(-0.11676107521566807256885,  0.42613108333816673845007, -1.24709647099828679905897, -1.59535854817539046734964,
                                           0.09466707711203618813567,  0.46023802097176974568171,  0.49675403299826370595582, -0.38191954281375029012224,
                                           -0.62952239136842547129191, -0.02095565193527237912718),
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
                  constraints.lambda = list(cbind(c(1,1,1,1))),
                  cor_struct_gauss = cor_ident())
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)
mixoglmm:::check(all.equal(params[, 1], c(1.565550788263897208363e-01, -1.005055677165627825431e+00, -1.155608614217086375919e-02,  9.067579782146614197913e-01,
                                          -9.449493971975334583036e-01, 4.674549843605066978824e-01,  1.000000000000000000000e+00,  6.642033527501332024201e-01,
                                          1.843609848911476311586e-05,  5.629183885736568360159e-01,  1.785009068471818283186e+00),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(params[, 2], c(0.07292725470527808284693, 0.06646600503762277456499, 0.07107379297895902170445, 0.10387610209229070079573,
                                          0.03229450514040876213384, 0.02997892096909285722384, 0.00000000000000000000000, 0.01395089666318996210470,
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
                cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_logit, method = "conditional means")

p_mode <- extract_ranef(fit_logit, method = "conditional modes")

ss <- summary(fit_logit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(-0.9546612843688029315103, -0.8482904907878927591725,  1.0000000000000000000000,
                                          0.1610632373419722818131, 2.8737288498619264487388,  0.6680940940787912785126),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.03739356462311806822418, 0.04399576042123849034526, 0.00000000000000000000000,
                                          0.04391550808809003419997, 1.74998817233197812015533, 0.03734077630716118451071),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(0.9309770137707257253723, 2.8064760580986050086949, -3.6623351637230201482964,
                                          -4.5272695540600462749126,-1.8731262869888314348543),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(0.877830973109587797687,  3.149437301772087582918, -3.852426316300728625919,
                                          -4.740299800718964284840, -2.202529149853830858774),
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
                      cor_struct_gauss = cor_ident(~ 1))

p_mean <- extract_ranef(fit_cloglog, method = "conditional means")

p_mode <- extract_ranef(fit_cloglog, method = "conditional modes")

ss <- summary(fit_cloglog)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(-0.9142100217791669924949, -0.9036882662885042938328,  1.0000000000000000000000,
                                          0.4196134074293071392070, 1.2937187108498517940092,  0.6008667556347018212648 ),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.03825668035180943699070, 0.04074373631613313617716, 0.00000000000000000000000,
                                          0.11785431500539235438740, 0.30917491040087202724962, 0.06346580389030342939094 ),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-0.3643462418018307835688,  1.2993685252252837880604, -2.0188697112039601933020,
                                          -2.5928120395154437005658, -1.0012653529242976357949),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-0.3978179417637880410652,  1.3654841614455761611424, -2.0445665210221735108576,
                                          -2.6342688498182487144561, -1.0480157963259864040140),
                           tol = tol, check.attributes = FALSE))


## log
# fit_log <- mixoglmm(formula2,
#                     families = list(
#                       Be1 = binomial(link="log"),
#                       y1 = gaussian()),
#                     data = data_toy,
#                     constraints.beta = list("(Intercept)" = cbind(c(1,1)),
#                                             X1 = cbind(c(1,1))),
#                     cor_struct_gauss = cor_general(~ 1))
#
# p_mean <- extract_ranef(fit_log, method = "conditional means")
#
# p_mode <- extract_ranef(fit_log, method = "conditional modes")

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
                        cor_struct_gauss = cor_general(~ 1))

p_mean <- extract_ranef(fit_cauchit, method = "conditional means")

p_mode <- extract_ranef(fit_cauchit, method = "conditional modes")

ss <- summary(fit_cauchit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.1451748792263497234156, 0.9491844792366731509148, -1.2421507284665880277430,
                                          0.5791218013574510603547, 1.0000000000000000000000,  0.2981999297653096792082,
                                          2.0072152644343770333535),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.10231743234024290312156, 0.13513258435161068504726, 0.08618636650438254998008,
                                          0.07155235029109439925943, 0.00000000000000000000000, 0.09530371279236171666582,
                                          0.08151635280000514571785),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mode[1:5], c(-0.00795667466004393482415,  0.04073846973298628909577, -0.07919674966990808417933,
                                          -0.13204059756172059625356, -0.04151928904885714838313),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(p_mean[1:5], c(-0.007806203806778966823587,  0.041045447858503336657776,
                                          -0.079989360189275729440261, -0.131372081678998414711046,
                                          -0.042438053741471443158773 ),
                           tol = tol, check.attributes = FALSE))



# ## Missing values
# withNA <- FALSE
# #data_toy$D1 <- as.numeric(data_toy$D1)-1
# #data_toy$firm <- seq_len(n)
# if (withNA) {
#   N <- nrow(data_toy)
#   set.seed(1000)
#   data_toy$Be1[sample(1:N, round(N/10))] <- NA
#   data_toy$Po1[sample(1:N, round(N/10))] <- NA
# }
#
#
# #summaryRprof("Rprof.out")
# object<-fit
# fit$y2
# fit
# print.summary.mixoglmm(summary(fit))
# vcov(fit)
#
# printCoefmat(summary(fit)$coefficients)
#
#
#
# cbind(pm1, pm2)
# fit$dims$N
#
#
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

### NAs
set.seed(12345)
data_toy[sample(1:nrow(data_toy), 50), "y1"] <- NA
data_toy[sample(1:nrow(data_toy), 50), "y2"] <- NA

## cor_general
formula <- (Be1 +  y1 + y2 + y3 ~ 1 + X1 + X2)
system.time(
  fit <- mixoglmm(formula, na.action = "na.pass",
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
                  cor_struct_gauss = cor_general(~ 1))
)
p_mean <- extract_ranef(fit, method = "conditional means")

p_mode <- extract_ranef(fit, method = "conditional modes")

options(digits = 22)
ss <- summary(fit)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.142523785053917939613655, -0.987446384577521429903868,  0.009187106859129847563628,
                                           0.925973973876359490731147, -0.952424773550239689434704,  0.467081141234807317719202,
                                           1.000000000000000000000000,  0.484396075957217286944712,  0.893731621943977816435734,
                                           0.501131179116041658438974,  0.447111238608769390179276,  0.445861707713270050135179,
                                           0.929225110144207167017782,  1.956983243496172031328229),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.06957874648604583933675, 0.06966806120683147196537, 0.07941034885812886523482,
                                        0.11057065074712447982908, 0.02917660515428029469054, 0.02684484294836699327935,
                                        0.00000000000000000000000, 0.02176931580807565935753, 0.01961876601887017901316,
                                        0.05028622805157003522458, 0.04171194856407005324694, 0.05153823987805021544784,
                                        0.04379678299428022697493, 0.06568141855243250781804 ),
                           tol = tol, check.attributes = FALSE))
## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode[1:5], c(0.26889256869879635258869, 0.36040527552096185415209,  0.81848672972749270115855,
                                           -0.66584384659365170033851,  0.09105434453282762463644),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean[1:5], c(0.26921503382808120719005, 0.36088204624612002513473,
                                          0.81810677143196663951841, -0.66629162975553590797517,  0.09035694045716302635896 ),
                           tol = tol, check.attributes = FALSE))
