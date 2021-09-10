#  
#     Script: Meta-analysis using metafor package
#     Aim: estimate response of pests (intrinsic rate of increase) to temperature
#          across categorical variables (taxa, feeding guilds, host plant, etc)


#### Prev. LOAD ####
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(metafor)
#### 1. Example Hamann 2021 ####
herbivore_data<-read.delim("C:/Users/Ecologia/Desktop/DARÍO_actualizada_30_06_2021/Intrinsic_metaanalysis/herbivore_data.txt", header=TRUE)
colnames(herbivore_data)
temp_hamann <- herbivore_data %>%
  filter(treatment=="Temperature") %>%
  select(study,trait,treatment,mean_control,SE_control,n_control,mean_treatment,SE_treatment,n_treatment)



#### 2. Our data ####
IR_data_metafor_raw<-read_csv("/Users/Ecologia/Desktop/DARÍO_actualizada_30_06_2021/Intrinsic_metaanalysis/IR_metaanalysis_suitable_arranged4metafor.csv")
acari <- IR_data_metafor_raw %>%
  filter(order == "Acari>Prostigmata" | 
           order == "Acari>Trombidiformes")%>%
  mutate(order = "Acari")
IR_data_metafor_raw2 <- IR_data_metafor_raw %>%
  filter(order != "Acari>Prostigmata" &
           order != "Acari>Trombidiformes")%>%
  bind_rows(acari)
# insert id
IR_id <- IR_data_metafor_raw2 %>%
  distinct(`Article Title`)%>%
  mutate(id=row_number())
IR_data_metafor <- inner_join(IR_data_metafor_raw2,IR_id,by="Article Title")
#select interest variables
eff_size_prep <- IR_data_metafor %>%
  select(id,10:17,order,feeding_guild,lon,lat,Authors,Year)%>%
  mutate(sd1=treatment_se*sqrt(treatment_n))%>%
  mutate(sd2=control_se*sqrt(control_n))%>%
  mutate(Authors= word(Authors,1,sep = ","))
 
eff_sizes_extest<-eff_size_prep[1:25,]
# compute effect sizes
eff_sizes <- escalc(measure="SMD", 
                    m2i=control_mean,
                    m1i=treatment_mean, 
                    sd2i=sd2, 
                    sd1i=sd1, 
                    n2i=control_n, 
                    n1i=treatment_n, 
                    data=eff_sizes_extest)

#un randomEffects sin agrupar por estudio (mal hecho)
res <- rma.uni(yi, vi, data=eff_sizes, method="REML")

#forest plot
forest(res, xlim=c(-0.0001, 0.00001))

## mixed effects ##
# compute effect sizes 
eff_sizes <- escalc(measure="SMD", 
                    m2i=control_mean,
                    m1i=treatment_mean, 
                    sd2i=sd2, 
                    sd1i=sd1, 
                    n2i=control_n, 
                    n1i=treatment_n, 
                    data=eff_size_prep)

eff_sizes_filtered <- eff_sizes %>%
  filter(yi<10 & vi <50)

res2 <- rma(yi, vi, data=eff_sizes,mods = id)
predict(res,newmods=cbind(seq(10,60,10)))
#forest plot
forest(res2, xlim=c(-0.0001, 0.00001))

## funnel plot
funnel(res,main="Mixed-Effects Model")

### subgroups ###
eff_sizes <- escalc(measure="SMD", 
                    m2i=control_mean,
                    m1i=treatment_mean, 
                    sd2i=sd2, 
                    sd1i=sd1, 
                    n2i=control_n, 
                    n1i=treatment_n, 
                    data=eff_size_prep)
res_lepi <- rma(yi, vi, subset=(order=="Lepidoptera"), data=eff_sizes)
forest(res_lepi)

### subgroups no independent ###
res_studies <- rma.mv(yi,vi,random = ~1|id,data=eff_sizes,
                      subset = (order=="Lepidoptera"))
forest(res_studies)
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects models to some subsets
res_lepi <- rma(yi, vi, subset=(order=="Lepidoptera"), data=eff_sizes)
res_hemip <- rma(yi, vi, subset=(order=="Hemiptera"), data=eff_sizes)
res_diptera <- rma(yi, vi, subset=(order=="Diptera"), data=eff_sizes)
res_acari <- rma(yi, vi, subset=(order=="Acari"), data=eff_sizes)

### collect model estimates and corresponding variances
estimates <- c(coef(res_lepi), coef(res_hemip), coef(res_diptera), coef(res_acari))
variances <- c(vcov(res_lepi), vcov(res_hemip), vcov(res_diptera),vcov(res_acari))

### create vector with labels
labels <- c("Lepidoptera", "Hemiptera", "Diptera","Acari")

### forest plot groups
forest(estimates, variances, slab=labels)

###similar, but with study variation
### fit random-effects models to some subsets
res_lepi <- rma.mv(yi, vi, subset=(order=="Lepidoptera"),
                   random=~1|id,data=eff_sizes)
res_hemip <- rma.mv(yi, vi, subset=(order=="Hemiptera"),
                 random=~1|id,data=eff_sizes)
res_diptera <- rma.mv(yi, vi, subset=(order=="Diptera"),
                   random=~1|id,data=eff_sizes)
res_acari <- rma.mv(yi, vi, subset=(order=="Acari"),
                 random=~1|id,data=eff_sizes)

### collect model estimates and corresponding variances
estimates <- c(coef(res_lepi), coef(res_hemip), coef(res_diptera), coef(res_acari))
variances <- c(vcov(res_lepi), vcov(res_hemip), vcov(res_diptera),vcov(res_acari))

### create vector with labels
labels <- c("Lepidoptera", "Hemiptera", "Diptera","Acari")

### forest plot groups
forest(estimates, variances, slab=labels)

###similar, but with study variation for FEEDING GUILD
### fit random-effects models to some subsets
res_sucker <- rma.mv(yi, vi, subset=(feeding_guild=="sucker"),
                   random=~1|id,data=eff_sizes)
res_borer <- rma.mv(yi, vi, subset=(feeding_guild=="borer"),
                    random=~1|id,data=eff_sizes)
res_chewer <- rma.mv(yi, vi, subset=(feeding_guild=="chewer"),
                      random=~1|id,data=eff_sizes)
res_miner <- rma.mv(yi, vi, subset=(feeding_guild=="miner" |
                                      feeding_guild =="leafminer"),
                    random=~1|id,data=eff_sizes)

### collect model estimates and corresponding variances
estimates <- c(coef(res_sucker), coef(res_borer), coef(res_chewer), coef(res_miner))
variances <- c(vcov(res_sucker), vcov(res_borer), vcov(res_chewer),vcov(res_miner))

### create vector with labels
labels <- c("Suckers", "Borers", "Chewers","Miners")

### forest plot groups
forest(estimates, variances, slab=labels)

#### aggregate_by_study ####
### fit random-effects models to some subsets
# (by study?)
eff_sizes_arranged <- eff_sizes %>%
  arrange(order)
agg<-aggregate(eff_sizes_arranged, cluster=id,
          struct="ID")
res_agg <- rma.mv(yi, vi,random=~1|id,data=agg)

forest(res_agg)
## order
res_lepi <- rma.mv(yi, vi, subset=(order=="Lepidoptera"),
                   random=~1|id,data=agg)
res_hemip <- rma.mv(yi, vi, subset=(order=="Hemiptera"),
                    random=~1|id,data=agg)
res_diptera <- rma.mv(yi, vi, subset=(order=="Diptera"),
                      random=~1|id,data=agg)
res_acari <- rma.mv(yi, vi, subset=(order=="Acari"),
                    random=~1|id,data=agg)

### collect model estimates and corresponding variances
estimates <- c(coef(res_lepi), coef(res_hemip), coef(res_diptera), coef(res_acari))
variances <- c(vcov(res_lepi), vcov(res_hemip), vcov(res_diptera),vcov(res_acari))

### create vector with labels
labels <- c("Lepidoptera", "Hemiptera", "Diptera","Acari")

### forest plot groups
forest(estimates, variances, slab=labels)
# it remains the same...

#### let's use raw mean rather than SDM ####
eff_sizes <- escalc(measure="MD", 
                    m2i=control_mean,
                    m1i=treatment_mean, 
                    sd2i=sd2, 
                    sd1i=sd1, 
                    n2i=control_n, 
                    n1i=treatment_n, 
                    data=eff_size_prep)
agg<-aggregate(eff_sizes, cluster=id,
               struct="ID")
summary(eff_sizes$vi)
eff_sizes_clean <- agg%>%
  arrange(id)%>%
  filter(vi <0.01)
res_agg <- rma.mv(yi, vi,
                  random=~1|id,
                  mods=~Year,
                  data=eff_sizes_clean,
                  slab=paste("Study:",eff_sizes_clean$id,eff_sizes_clean$Authors,",",eff_sizes_clean$Year, sep=" "))
forest(res_agg,showweights = T,header = TRUE,)
### funnel plots ###
funnel(res_agg)
## order
res_lepi <- rma.mv(yi, vi, subset=(order=="Lepidoptera"),
                   random=~1|id,data=eff_sizes_clean)
res_hemip <- rma.mv(yi, vi, subset=(order=="Hemiptera"),
                    random=~1|id,data=eff_sizes_clean)
res_diptera <- rma.mv(yi, vi, subset=(order=="Diptera"),
                      random=~1|id,data=eff_sizes_clean)
res_acari <- rma.mv(yi, vi, subset=(order=="Acari"),
                    random=~1|id,data=eff_sizes_clean)

### collect model estimates and corresponding variances
estimates <- c(coef(res_lepi), coef(res_hemip), coef(res_diptera), coef(res_acari))
variances <- c(vcov(res_lepi), vcov(res_hemip), vcov(res_diptera),vcov(res_acari))

### create vector with labels
labels <- c("Lepidoptera", "Hemiptera", "Diptera","Acari")

### forest plot groups
forest(estimates, variances, slab=labels,showweights = T,header = TRUE)


#### Other APPROACH ###
#### _ _ 1. Using coefficients rather than mean diff ####
eff_sizes <- escalc(measure="ZCOR", 
                    m2i=control_mean,
                    m1i=treatment_mean, 
                    sd2i=sd2, 
                    sd1i=sd1, 
                    n2i=control_n, 
                    n1i=treatment_n, 
                    data=eff_size_prep)
