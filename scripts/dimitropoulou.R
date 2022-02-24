# example simulations for meta-analysis  --------- 
# from Papadimitropoulou et al. (2019) 
#(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6767371)

data1 <- data.frame(study = c(1:5,1:5),
                    mean = c(114,89,63,101,81,100,90,56,114,67),
                    sd = c(25,32,30,31,31,39,33,22,35,23),
                    n = c(26,20,421,28,50,20,26,31,26,50),
                    group = c(rep(0,5),rep(1,5)))
# generate IPD
# First expand dataset and generate obesrvations from a standard N(0,1) distribution
data_IPD <- data.frame(study_i = rep(data1$study, data1$n),
                       mean_i = rep(data1$mean, data1$n),
                       sd_i = rep(data1$sd, data1$n),
                       group_i = rep(data1$group, data1$n))
set.seed(64870236)
data_IPD$ytmp = rnorm(n = nrow(data_IPD),
                      0,
                      1)
# Now generate from Y*IPD data with exactly the same mean and sd as in the summary data
mean_IPDtmp <- aggregate(ytmp ~ group_i + study_i,
                         data_IPD, 
                         mean)
sd_IPDtmp <- aggregate(ytmp ~ group_i + study_i,
                         data_IPD, 
                         sd)
names(mean_IPDtmp)[names(mean_IPDtmp) == "ytmp"] = "ytmpmean"
names(sd_IPDtmp)[names(sd_IPDtmp) == "ytmp"] = "ytmpsd"

data_IPD = merge(data_IPD, mean_IPDtmp, by = c("study_i", "group_i"))
data_IPD = merge(data_IPD, sd_IPDtmp, by = c("study_i", "group_i"))
data_IPD$y = data_IPD$mean_i + (data_IPD$ytmp - data_IPD$ytmpmean)*(data_IPD$sd_i/data_IPD$ytmpsd)
data_IPD$arm = 1000 * data_IPD$study_i + data_IPD$group_i
view(data_IPD)
# with dplyr
sum_dplyr<- data_IPD %>% 
  group_by(group_i, study_i) %>% 
  summarise(ytmpmean = mean(ytmp),
            ytmpsd = sd(ytmp)) 

merging <- inner_join(data_IPD,sum_dplyr) %>% 
  mutate(y = mean_i + (ytmp - ytmpmean)*(sd_i/ytmpsd),
         arm = 1000 * study_i + group_i)  
write.table(data1, file = "dimitropoulou_raw.txt", sep = ",", quote = FALSE, row.names = F)
write.table(merging, file = "dimitropoulou.txt", sep = ",", quote = FALSE, row.names = F)
