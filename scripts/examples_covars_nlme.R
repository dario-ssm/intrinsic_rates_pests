#### examples nlme covariates ####
library(nlme)
## co2
CO2
fm2CO2 <- nlme(uptake ~ SSasympOff(conc, Asym, lrc, c0),
               data = CO2,
               fixed = Asym + lrc+ c0 ~ 1,
               random = Asym + lrc ~ 1)
fm2CO2
ranef(fm2CO2, augFrame = T)
fm3C02 <- update(fm2CO2,
                 fixed = list(Asym + lrc ~ Type * Treatment, c0 ~ 1),
                 start = c(32.4, 0, 0, 0, -4.6, 0, 0, 0, 49.3))
fm3C02

## myData
starts_all_sim_sens <-c(a = 9e-05, Tmin = 9, Tmax = 40)
IR_sim_ord_subset_last <- IR_sim_ord_subset %>% 
  mutate(taxa = as.factor(order)) %>% 
  filter(id != 37)
nlme_br1_order_subset_test <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                                   fixed = a+Tmin+Tmax ~ 1,
                                   start = starts_all_sim_sens,
                                   groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                                   data = IR_sim_ord_subset_last,
                                   weights = varExp(),
                                   control = nlmeControl(pnlsTol = 1,
                                                         msMaxIter = 100,
                                                         maxIter = 100,
                                                         msVerbose = TRUE))
sum_nlme <- summary(nlme_br1_order_subset_test)
starts_nlme <- sum_nlme$tTable[,1]
nlme_br1_order_subset_test <- nlme(model= r ~ briere1(a,temp,Tmin,Tmax),
                              fixed = a+Tmin+Tmax ~ taxa,
                              start = c(starts_nlme[1], 0, 0, # <- a
                                        starts_nlme[2], 0, 0, # <- Tmin,
                                        starts_nlme[3], 0, 0 # <- Tmax
                                        ),
                              groups = ~id, # same as random effects = a+Tmin+Tmax ~ 1| id
                              data = IR_sim_ord_subset_last,
                              weights = varExp(),
                              control = nlmeControl(pnlsTol = 10,
                                                    msMaxIter = 50,
                                                    maxIter = 50,
                                                    msVerbose = TRUE))
summary(nlme_br1_order_subset_test)
            



   ### bayes??

prior_1 <- c(set_prior("normal(1.317e-04, (2.855e-05)^2)", nlpar = "a"),
             set_prior("normal(8.158e+00, (2.546e+00)^2)", nlpar = "Tmin"),
             set_prior("normal(4.097e+01, (1.979e+00)^2)", nlpar = "Tmax"))

briere_formula <- brmsformula(formula =  r ~ a*temp*(temp-Tmin)*(Tmax-temp)^(1/2),
                              # Nonlinear variables
                              nonlinear = list(a+Tmin+Tmax ~ 1|id),
                              # Nonlinear fit
                              nl = TRUE)
bayes_fit <- brm(
  briere_formula,
  family=gaussian(), 
  data = acari_data,
  prior = prior_1)
summary(bayes_fit)


n2_b <- brm(bf(f1, nonlinear = list(phi1 ~ (1|Tree),
                                    phi2 ~ 1,
                                    phi3 ~ 1),nl=TRUE),
            data = Orange,
            prior = prior_2,
            chains = 3,
            iter = 4000)