#varIdent por id
ident_nls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = c(a = 9e-05, Tmin = 8, Tmax = 40),
                  weights = varIdent(form = ~1|id), 
                  data = intrapest_test)
sum_ident <- summary(ident_nls)
ident_vals <- c(sum_ident$coefficients, sum_ident$logLik, sum_ident$AIC)
#ident por vi
ident_nls2 <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                   start = c(a = 9e-05, Tmin = 8, Tmax = 40),
                   weights = varIdent(form = ~sigma), 
                   data = intrapest_test)
sum_ident2 <- summary(ident_nls2)
ident_vals2 <- c(sum_ident2$coefficients, sum_ident2$logLik, sum_ident2$AIC)


#varExp
exp_nls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = c(a = 9e-05, Tmin = 8, Tmax = 40),
                  weights = varExp(form = ~temp), 
                  data = intrapest_test)
sum_exp <- summary(exp_nls)
exp_vals <- c(sum_exp$coefficients, sum_ident$logLik, sum_ident$AIC)

#comb(varIdent por id, varExp)
comb_nls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                start = c(a = 9e-05, Tmin = 8, Tmax = 40),
                weights = varComb(varIdent(form = ~ 1|id),
                                  varExp(form = ~ temp)),
                data = intrapest_test,
                control = gnlsControl(nlsTol = 1, nlsMaxIter = 15))
sum_comb <- summary(comb_nls)
comb_vals <- c(sum_comb$coefficients, sum_comb$logLik, sum_comb$AIC)

#comb(varIdent por vi, varExp)
comb2_gnls <- gnls(r ~ briere1(a,temp = temp,Tmin,Tmax),
                  start = c(a = 9e-05, Tmin = 8, Tmax = 40),
                  weights = varComb(varIdent(form = ~ sigma),
                                    varExp(form = ~ temp)),
                  data = intrapest_test,
                  control = gnlsControl(nlsTol = 1, nlsMaxIter = 15))
sum_comb2 <- summary(comb2_gnls)
comb_vals2 <- c(sum_comb2$coefficients, sum_comb2$logLik, sum_comb2$AIC)






names_comparative <- c("a","Tmin","Tmax","log-lik","AIC")
comparative <- rbind(ident_vals, ident_vals2, exp_vals, comb_vals, comb_vals2)
colnames(comparative) <- names_comparative
view(comparative)

v1 <- c("A", "B", "C")
vrep <- rep(v1, each = 3)
