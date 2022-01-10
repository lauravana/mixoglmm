tol <- 1e-4
####################
devtools::install_github("lauravana/mixoglmm")
library(mixoglmm)
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

# paste0(params[, 1], collapse = ", ")
mixoglmm:::check(all.equal(params[, 1], c(-0.954662618387938, -0.848290475403047, 1,
                                          0.184391890658427, 2.97914016287556, 0.714841872699901),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.037393571603032, 0.0439957011026461, 0, 0.0479925611997275,
                                          1.85470557715777, 0.0431556051459695),
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

mixoglmm:::check(all.equal(params[, 1], c(-0.914210189304185, -0.903688665092929, 1
                                          , 0.46870167199286, 1.36029574788433, 0.635041794662754),
                        tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.0382566560100977, 0.0407437133603616, 0,
                                          0.12137012853852, 0.328398892202561, 0.0740727803639353),
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

mixoglmm:::check(all.equal(params[, 1], c(0.149226223992418721709, 0.945156020213949732423, -1.241323798981201598224,
                                          0.579684298052147384261, 1.0000000000000000000000,  0.321700803542295687976,
                                          1.388005722277994236080),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[, 2], c(0.101942986184902661950, 0.133205452009057229157, 0.085435873751396690379,
                                          0.072062241895121767477, 0.00000000000000000000000, 0.051612301832824669656,
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
                  cor_struct_gauss = cor_general(~ 1))
)

p_mean_NA <- extract_ranef(fit_NA, method = "conditional means")

p_mode_NA <- extract_ranef(fit_NA, method = "conditional modes")

options(digits = 22)
ss <- summary(fit_NA)
params <- rbind(ss$coefficients, ss$re.coefficients, ss$re.stddev, ss$gauss.corr, ss$gauss.stddev)

mixoglmm:::check(all.equal(params[, 1], c(0.153623278074614, -0.99951275409134, -0.0145107052108456, 0.910636241262308, -0.952695650977913, 0.475271346405206,
                                          1, 0.459786920277783, 0.89776837213671, 0.50175237066725, 0.453896609941155,
                                          0.472387980292025, 0.954880636610505, 1.97329841750088),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(params[,2], c(0.0694463900632229, 0.0700169815010208, 0.0798689916318538,
                                         0.110858992073995, 0.0287371190721126, 0.0265436389443322,
                                         0, 0.0205925743953325, 0.0165299203289168,
                                         0.0476127317388988, 0.0408248560466499, 0.0478621465305282,
                                         0.0436891254674255, 0.0660832026795537),
                           tol = tol, check.attributes = FALSE))

tau <- ss$re.stddev
mixoglmm:::check(all.equal(tau[, 1], c(0.4597869202777826180828),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(tau[, 2], c(0.02059257439533253644659 ),
                           tol = tol, check.attributes = FALSE))

gamma <- ss$gauss.corr
mixoglmm:::check(all.equal(gamma[, 1], c(0.897768372136710457454, 0.5017523706672499805848, 0.4538966099411552712617),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(gamma[, 2], c(0.01652992032891676885131, 0.04761273173889878779219, 0.04082485604664993356083),
                           tol = tol, check.attributes = FALSE))

omega <- ss$gauss.stddev
mixoglmm:::check(all.equal(omega[,1], c(0.4723879802920253889731, 0.9548806366105048182291, 1.9732984175008785321381),
                           tol = tol, check.attributes = FALSE))
mixoglmm:::check(all.equal(omega[,2], c(0.04786214653052819900658, 0.04368912546742549640744, 0.06608320267955368476631),
                           tol = tol, check.attributes = FALSE))


## checked by hand with analytical for model, OK
mixoglmm:::check(all.equal(p_mode_NA[1:5], c(0.0168590889448626, 0.328255749701425, -0.554991843035224, -0.595085120283989, 0.0992605827140897),
                           tol = tol, check.attributes = FALSE))

mixoglmm:::check(all.equal(p_mean_NA[1:5], c(0.0171581467395677, 0.328798850431899, -0.555755506012409, -0.595577552428667, 0.0984938254729672),
                           tol = tol, check.attributes = FALSE))
