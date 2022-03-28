# SCRIPT INFO --------------------
#     Authors: Dario San Segundo Molina, Sara Villen Perez, Ignacio Morales Castilla
#     Title: meta-analyses plots
#     Aim: plot results of meta-analytical models for the dataset of thermal traits
#     Date: March 2022
#
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
library(tidyverse)
library(ggthemes)
library(magrittr)
library(wesanderson)
library(viridis)
library(hrbrthemes)
library(networkD3)

# 1. Calculate i-th study uncertainties ----
## note that no transformation is required here since these intervals are being used
## as an approximation for forest plotting
### complete intervals
intervals_i <- thermal_traits_complete %>% 
  mutate(tmax_ci_lower = tmax-1.96*tmax_se,
         tmax_ci_upper = tmax +1.96*tmax_se,
         tmin_ci_lower = tmin-1.96*tmin_se,
         tmin_ci_upper = tmin +1.96*tmin_se,
         topt_ci_lower = topt-1.96*topt_se,
         topt_ci_upper = topt +1.96*topt_se) %>% 
  select(id, tmax, tmax_se, tmin, tmin_se, topt, topt_se, thermal_breadth,vi,
         tmax_ci_lower,
         tmax_ci_upper,
         tmin_ci_lower,
         tmin_ci_upper,
         topt_ci_lower,
         topt_ci_upper) %>% 
  group_by(id) %>% 
  summarise_all(mean) %>% 
  print()

### include ranges
ranges_i <- thermal_traits_complete %>% 
  group_by(id) %>% 
  summarise(tmax_range_lower = range(tmax)[1],
            tmax_range_upper = range(tmax)[2],
            tmin_range_lower = range(tmin)[1],
            tmin_range_upper = range(tmin)[2],
            topt_range_lower = range(topt)[1],
            topt_range_upper = range(topt)[2],
            thermal_breadth_range_lower = range(thermal_breadth)[1],
            thermal_breadth_range_upper = range(thermal_breadth)[2]) %>% 
  select(-id)
authors_i <- thermal_traits_complete %>% 
  group_by(id) %>% 
  summarise(authors = unique(authors),
            year = mean(Year),
            order = unique(order)) %>%
  mutate(authors = paste0(authors,", ",year)) %>% 
  select(-id)
uncertainty_i <- intervals_i %>% 
  bind_cols(ranges_i,
            authors_i) %>% 
  select(-year) %>% 
  mutate(ggweights = 1/vi)

# 1. Forest plots summary  ----
# ~~ a) tmax  ----
## add summary effect and uncertainty
sum_eff_tmax_int <- as_tibble(tmax_intercept_output) 
forest_tmax_sum <- ggplot(data = uncertainty_i, aes(x = tmax, y = fct_reorder(as_factor(id),order)))+
  geom_pointrange(aes(xmin = tmax_ci_lower,
                      xmax = tmax_ci_upper,
                      color = order))+
  geom_point(aes(size = ggweights, color = order))+
  scale_y_discrete(expand = c(0,2)) +
  scale_color_brewer(palette = "Set2")+
  theme_few()+
  annotate(geom = "segment",
           x = sum_eff_tmax_int$ci_l,
           xend = sum_eff_tmax_int$ci_h,
           y = 0,
           yend = 0,
           size = 1.5,
           color = "deeppink4")+
  annotate(geom = "rect", xmin = 34.1, xmax = 34.9,
           ymin = -0.5, ymax = 0.5, fill = "deeppink4")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
        )+
  labs(x = "Maximum Temperature (ºC)",
       y = "")
forest_tmax_sum
ggsave(filename = "forest_tmax_sumeff.png", width = 20, height = 30, dpi = 300, units = "cm")
# ~~ b) tmin  ----
## add summary effect and uncertainty
sum_eff_tmin_int <- as_tibble(tmin_intercept_output) 
forest_tmin_sum <- ggplot(data = uncertainty_i, aes(x = tmin, y = fct_reorder(as_factor(id),order)))+
  geom_pointrange(aes(xmin = tmin_ci_lower,
                      xmax = tmin_ci_upper,
                      color = order))+
  geom_point(aes(size = ggweights, color = order))+
  scale_y_discrete(expand = c(0,3)) +
  scale_color_brewer(palette = "Set2")+
  theme_few()+
  annotate(geom = "segment",
           x = sum_eff_tmin_int$ci_l,
           xend = sum_eff_tmin_int$ci_h,
           y = -0.5,
           yend = -0.5,
           size = 1.5,
           color = "deeppink4")+
  annotate(geom = "rect", xmin = sum_eff_tmin_int$pred-0.3, xmax = sum_eff_tmin_int$pred+0.3,
           ymin = -1, ymax = 0, fill = "deeppink4")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )+
  labs(x = "Minimum Temperature (ºC)",
       y = "")
forest_tmin_sum
ggsave(filename = "forest_tmin_sumeff.png", width = 20, height = 30, dpi = 300, units = "cm")

# ~~ c) topt  ----
## add summary effect and uncertainty
sum_eff_topt_int <- as_tibble(topt_intercept_output) 
forest_topt_sum <- ggplot(data = uncertainty_i, aes(x = topt, y = fct_reorder(as_factor(id),order)))+
  geom_pointrange(aes(xmin = topt_ci_lower,
                      xmax = topt_ci_upper,
                      color = order))+
  geom_point(aes(size = ggweights, color = order))+
  scale_y_discrete(expand = c(0,3)) +
  scale_color_brewer(palette = "Set2")+
  theme_few()+
  annotate(geom = "segment",
           x = sum_eff_topt_int$ci_l,
           xend = sum_eff_topt_int$ci_h,
           y = -0.5,
           yend = -0.5,
           size = 1.5,
           color = "deeppink4")+
  annotate(geom = "rect", xmin = sum_eff_topt_int$pred-0.3, xmax = sum_eff_topt_int$pred+0.3,
           ymin = -1, ymax = 0, fill = "deeppink4")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )+
  labs(x = "Optimal Temperature (ºC)",
       y = "")
forest_topt_sum
ggsave(filename = "forest_topt_sumeff.png", width = 20, height = 30, dpi = 300, units = "cm")

# ~~ d) thermal_breadth  ----
## add summary effect and uncertainty
sum_eff_thermal_breadth_int <- as_tibble(thermal_breadth_intercept_output) 
forest_thermal_breadth_sum <- ggplot(data = uncertainty_i, aes(x = thermal_breadth, y = fct_reorder(as_factor(id),order)))+
  geom_pointrange(aes(xmin = thermal_breadth_ci_lower,
                      xmax = thermal_breadth_ci_upper,
                      color = order))+
  geom_point(aes(size = ggweights, color = order))+
  scale_y_discrete(expand = c(0,3)) +
  scale_color_brewer(palette = "Set2")+
  theme_few()+
  annotate(geom = "segment",
           x = sum_eff_thermal_breadth_int$ci_l,
           xend = sum_eff_thermal_breadth_int$ci_h,
           y = -0.5,
           yend = -0.5,
           size = 1.5,
           color = "deeppink4")+
  annotate(geom = "rect", xmin = sum_eff_thermal_breadth_int$pred-0.3, xmax = sum_eff_thermal_breadth_int$pred+0.3,
           ymin = -1, ymax = 0, fill = "deeppink4")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )+
  labs(x = "Thermal Breadth (?C)",
       y = "")
forest_thermal_breadth_sum
# ~~ d) thermal safety margin  ----
## add summary effect and uncertainty
sum_eff_tsm_int <- as_tibble(tsm_intercept_output) 
forest_tsm_sum <- ggplot(data = uncertainty_i, aes(x = tsm, y = fct_reorder(as_factor(id),order)))+
  geom_pointrange(aes(xmin = tsm_ci_lower,
                      xmax = tsm_ci_upper,
                      color = order))+
  geom_point(aes(size = ggweights, color = order))+
  scale_y_discrete(expand = c(0,3)) +
  scale_color_brewer(palette = "Set2")+
  theme_few()+
  annotate(geom = "segment",
           x = sum_eff_tsm_int$ci_l,
           xend = sum_eff_tsm_int$ci_h,
           y = -0.5,
           yend = -0.5,
           size = 1.5,
           color = "deeppink4")+
  annotate(geom = "rect", xmin = sum_eff_tsm_int$pred-0.3, xmax = sum_eff_tsm_int$pred+0.3,
           ymin = -1, ymax = 0, fill = "deeppink4")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )+
  labs(x = "Thermal Safety Margin (?C)",
       y = "")
forest_tsm_sum

# 2. Forest plots order  ----
# ~~ a) thermal points  ----
## add summary effect and uncertainty
sum_eff_tmax_order <- as_tibble(tmax_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "tmax")
sum_eff_tmin_order <- as_tibble(tmin_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "tmin")
sum_eff_topt_order <- as_tibble(topt_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "topt")

wes1 <- wes_palettes$Zissou1

traits_sumeff_order <- sum_eff_tmax_order %>% 
  bind_rows(sum_eff_tmin_order,
            sum_eff_topt_order)
forest_traits_order <- ggplot(data = traits_sumeff_order, aes(x = pred, y = order))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = trait),
                  size = 1.12,
                  position = position_jitter(height = 0.02))+
  scale_color_manual(values = wes1[c(5,1,3)])+
  theme_few()+
  labs(x = "Temperature (ºC)",
       y = "Order")
forest_traits_order
ggsave(filename = "forest_traits_order.png", width = 20, height = 16, dpi = 300, units = "cm")

# ~~ b) thermal_breadth  ----
## add summary effect and uncertainty
wes2 <- wes_palettes$Darjeeling1
sum_eff_thermal_breadth_order <- as_tibble(thermal_breadth_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "thermal breadth")
forest_thermal_breadth_order <- ggplot(data = sum_eff_thermal_breadth_order, 
                                       aes(x = pred, y = order))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = order),
                  size = 1.12)+
  scale_color_manual(values = wes2[c(1,2,3)])+
  theme_few()+
  labs(x = "Thermal Breadth (Δ ºC)",
       y = "Order")
forest_thermal_breadth_order
ggsave(filename = "forest_thermal_breadth_order.png",  width = 20, height = 16, dpi = 300, units = "cm")

# ~~ c) thermal safety margin  ----
## add summary effect and uncertainty
wes3 <- wes_palettes$BottleRocket2
sum_eff_tsm_order <- as_tibble(tsm_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "thermal breadth")
forest_tsm_order <- ggplot(data = sum_eff_tsm_order,
                           aes(x = pred, y = order))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = order),
                  size = 1.12)+
  scale_color_manual(values = wes3[c(1,2,3)])+
  theme_few()+
  labs(x = "Thermal Safety Margin (Δ ºC)",
       y = "Order")
forest_tsm_order
ggsave(filename = "forest_tsm_order.png",  width = 20, height = 16, dpi = 300, units = "cm")


# ~~ d) thermal range  ----
## add summary effect and uncertainty
wes4 <- wes_palettes$Cavalcanti1
sum_eff_therm_range_order <- as_tibble(therm_range_order_output) %>% 
  inner_join(counts_by_order, by = "order") %>% 
  mutate(trait = "thermal breadth")
forest_therm_range_order <- ggplot(data = sum_eff_therm_range_order,
                           aes(x = pred, y = order))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = order),
                  size = 1.12)+
  scale_color_manual(values = wes4[c(1,4,5)])+
  theme_few()+
  labs(x = "Thermal Range (Δ ºC)",
       y = "Order")
forest_therm_range_order
ggsave(filename = "forest_therm_range_order.png", width = 20, height = 16, dpi = 300, units = "cm")

# 3. Forest plots feeding_guild  ----
# ~~ a) thermal points  ----
## add summary effect and uncertainty
sum_eff_tmax_feeding_guild <- as_tibble(tmax_feeding_guild_output) %>% 
  mutate(trait = "tmax")
sum_eff_tmin_feeding_guild <- as_tibble(tmin_feeding_guild_output) %>% 
  mutate(trait = "tmin")
sum_eff_topt_feeding_guild <- as_tibble(topt_feeding_guild_output) %>% 
  mutate(trait = "topt")

wes1 <- wes_palettes$Zissou1

traits_sumeff_feeding_guild <- sum_eff_tmax_feeding_guild %>% 
  bind_rows(sum_eff_tmin_feeding_guild,
            sum_eff_topt_feeding_guild)
forest_traits_feeding_guild <- ggplot(data = traits_sumeff_feeding_guild, aes(x = pred, y = feeding_guild))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = trait),
                  size = 1.12,
                  position = position_jitter(height = 0.02))+
  scale_color_manual(values = wes1[c(5,1,3)])+
  theme_few()+
  labs(x = "Temperature (ºC)",
       y = "Feeding guild")
forest_traits_feeding_guild
ggsave(filename = "forest_traits_feeding_guild.png", width = 16, height = 16, dpi = 300, units = "cm")

# ~~ b) thermal_breadth  ----
## add summary effect and uncertainty
wes2 <- wes_palettes$Darjeeling1
sum_eff_thermal_breadth_feeding_guild <- as_tibble(thermal_breadth_feeding_guild_output)  %>% 
  mutate(trait = "thermal breadth")
forest_thermal_breadth_feeding_guild <- ggplot(data = sum_eff_thermal_breadth_feeding_guild, 
                                       aes(x = pred, y = feeding_guild))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = feeding_guild),
                  size = 1.12)+
  scale_color_manual(values = wes2[c(1,2,3)])+
  theme_few()+
  labs(x = "Thermal Breadth (Δ ºC)",
       y = "Feeding guild")+
  theme(legend.position = "none")
forest_thermal_breadth_feeding_guild
ggsave(filename = "forest_thermal_breadth_feeding_guild.png", width = 16, height = 16, dpi = 300, units = "cm")

# ~~ c) thermal safety margin  ----
## add summary effect and uncertainty
wes3 <- wes_palettes$BottleRocket2
sum_eff_tsm_feeding_guild <- as_tibble(tsm_feeding_guild_output) %>% 
  mutate(trait = "thermal breadth")
forest_tsm_feeding_guild <- ggplot(data = sum_eff_tsm_feeding_guild,
                           aes(x = pred, y = feeding_guild))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = feeding_guild),
                  size = 1.12)+
  scale_color_manual(values = wes3[c(1,2,3)])+
  theme_few()+
  labs(x = "Thermal Safety Margin (Δ ºC)",
       y = "Feeding guild")+
  theme(legend.position = "none")
forest_tsm_feeding_guild
ggsave(filename = "forest_tsm_feeding_guild.png", width = 16, height = 16, dpi = 300, units = "cm")

# ~~ d) thermal range  ----
## add summary effect and uncertainty
wes4 <- wes_palettes$Cavalcanti1
sum_eff_therm_range_feeding_guild <- as_tibble(therm_range_feeding_guild_output) %>% 
  mutate(trait = "thermal breadth")
forest_therm_range_feeding_guild <- ggplot(data = sum_eff_therm_range_feeding_guild,
                                   aes(x = pred, y = feeding_guild))+
  geom_pointrange(aes(xmin = ci_l,
                      xmax = ci_h,
                      color = feeding_guild),
                  size = 1.12)+
  scale_color_manual(values = wes4[c(1,4,5)])+
  theme_few()+
  labs(x = "Thermal Range (ΔºC)",
       y = "Feeding guild")+
  theme(legend.position = "none")
forest_therm_range_feeding_guild
ggsave(filename = "forest_therm_range_feeding_guild.png", width = 16, height = 16, dpi = 300, units = "cm")

# ~~ e) Relationship order - feeding guild ----
# Now we have 2 data frames: a 'links' data frame with 3 columns (from, to, value), and a 'nodes' data frame that gives the name of each node.
structures <- thermal_traits_complete %>%
  select(id, order, feeding_guild, genus, species) %>% 
  group_by(id) %>% 
  summarise_all(unique) %>% 
  mutate(spp = paste(genus, species)) %>% 
  select(order, spp)
simpleNetwork(structures, nodeColour = wes1)
# make a nodes data frame out of all unique nodes in networkData
nodes <- data.frame(name = unique(c(structures$order, structures$spp)))
# make a group variable where nodes in networkData$src are identified
nodes$group <- nodes$name %in% structures$order
# make a links data frame using the indexes (0-based) of nodes in 'nodes'
links <- data.frame(source = match(structures$order, nodes$name) - 1,
                    target = match(structures$spp, nodes$name) - 1)

forceNetwork(Links = links, 
             Nodes = nodes,
             Source = "source",
             Target = "target",
             NodeID ="name",
             Group = "group",
             zoom = 100,
             opacity = 1, opacityNoHover = 0.8)

# save the widget
# library(htmlwidgets)
# saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyEnergy.html"))
sankeyNetwork(Links = paramete, Nodes = Energy$nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30)


structures2 <- thermal_traits_complete %>%
  select(id, order, feeding_guild, genus, species) %>% 
  group_by(id) %>% 
  summarise_all(unique) %>% 
  mutate(spp = paste(genus, species)) %>% 
  select(order, feeding_guild) %>% 

simpleNetwork(structures2, nodeColour = wes1)
# make a nodes data frame out of all unique nodes in networkData
nodes <- data.frame(name = unique(c(structures2$feeding_guild, structures2$order)))
# make a group variable where nodes in networkData$src are identified
nodes$group <- nodes$name %in% structures2$feeding_guild
# make a links data frame using the indexes (0-based) of nodes in 'nodes'
links <- data.frame(source = match(structures2$feeding_guild, nodes$name) - 1,
                    target = match(structures2$order, nodes$name) - 1)
structures_both <- igraph_to_networkD3(structures2)
forceNetwork(Links = links, 
             Nodes = nodes,
             Source = "source",
             Target = "target",
             NodeID ="name",
             Group = "group",
             zoom = 100,
             opacity = 1, opacityNoHover = 0.8,
             fontSize = 14, bounded = TRUE) %>% 
  saveNetwork(file = 'Network1.html')

# 4. Feeding guild composition  ----
# ~~ a) barplots  ----
composition <- suckers_composition %>% 
  bind_rows(borers_composition, chewers_composition) %>% 
  print()
composition_bar <- ggplot(composition, aes(x = feeding_guild, y = perc, fill = order))+
                          geom_bar(position="stack", stat = "identity", color = "transparent")+
                          scale_fill_brewer(palette =  "Dark2") +
                          ggtitle("Feeding guild taxonomic composition")+
                          theme_fivethirtyeight()
  
ggsave("composition_bar.png", width = 30, height = 30, units = "cm")



# 5. Thermal traits across latitude  ----
## ~~ a) thermal points ----
Tmin_to_lat_plot <- ggplot(thermal_traits_complete,
                           aes(abs(lat),tmin))+
  geom_point(alpha=0.032,
             color="skyblue4")+
  geom_smooth(method="lm",color="skyblue4",fill="skyblue2")+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  theme_few()
tmin_to_lat_plot 

tmin_to_lat_order_plot <- ggplot(thermal_traits_order,aes(abs(lat),tmin))+
  geom_point(alpha=0.032,
             aes(color = order))+
  geom_smooth(method="lm", aes(color = order, fill=order))+
  scale_color_manual(values = wes2[c(2,3,5)])+
  scale_fill_manual(values = wes2[c(2,3,5)])+
  labs(title= "Tmin ~ latitude",
       x= "latitude",
       y= "Tmin")+
  facet_wrap(.~order)+
  theme_few()
tmin_to_lat_order_plot

## tmax
tmax_to_lat_plot <- ggplot(thermal_traits_complete,aes(abs(lat),tmax))+
  geom_point(alpha=0.032,
             color="red4")+
  geom_smooth(method="lm",color="red4",fill="red3")+
  labs(title= "Tmax ~ latitude",
       x= "latitude",
       y= "Tmax")+
  theme_classic()
tmax_to_lat_plot
tmax_to_lat_order_plot <- ggplot(thermal_traits_order,aes(abs(lat),tmax))+
  geom_point(alpha=0.032,
             aes(color = order))+
  geom_smooth(method="lm", aes(color = order, fill=order))+
  scale_color_manual(values = wes2[c(2,3,5)])+
  scale_fill_manual(values = wes2[c(2,3,5)])+
  labs(title= "tmax ~ latitude",
       x= "latitude",
       y= "tmax")+
  facet_wrap(.~order)+
  theme_few()
tmax_to_lat_order_plot

## topt
topt_to_lat_plot <- ggplot(thermal_traits_complete,aes(abs(lat),topt))+
  geom_point(alpha=0.032,
             color="mediumorchid4")+ 
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  theme_classic()
topt_to_lat_plot
topt_to_lat_order_plot <- ggplot(thermal_traits_complete,aes(abs(lat),topt))+
  geom_point(alpha=0.032,
             color="mediumorchid4")+
  geom_smooth(method="lm",color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Topt ~ latitude",
       x= "latitude",
       y= "Topt")+
  facet_wrap(.~order)+
  theme_light()
topt_to_lat_order_plot

lms_plot_grid <- cowplot::plot_grid(Tmin_to_lat_plot,Tmax_to_lat_plot,
                           nrow = 1,labels = c("A","B"))
lms_plot_grid
ggsave("tmin_tmax_lat_plotgrid.png", width = 24, height = 16, units = "cm", dpi = 300)

all_loess_combined <- ggplot(thermal_traits_complete)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3")+
  geom_point(aes(x=abs(lat),y=topt),color="mediumorchid4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=topt),color="mediumorchid4",fill="mediumorchid3")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()
all_loess_combined  

all_lms_combined <- ggplot(thermal_traits_complete)+
  geom_point(aes(x=abs(lat),y=tmin),color="skyblue4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmin),color="skyblue4",fill="skyblue2", method = "lm")+
  geom_point(aes(x=abs(lat),y=tmax),color="red4",alpha=0.032)+
  geom_smooth(aes(x=abs(lat),y=tmax),color="red4",fill="red3", method = "lm")+
  labs(title= "Thermal traits ~ latitude",x= "latitude",y= "Thermal traits (ºC)")+
  theme_classic()
all_lms_combined  
ggsave("tmin_tmax_lat_unique.png", width = 16, height = 16, units = "cm")

## ~~ b) thermal breadth ----
breadth_lat <- ggplot(data = thermal_traits_complete, aes(abs(lat), thermal_breadth))+
  geom_point(color = wes2[5], alpha = 0.1)+
  geom_smooth(method = "lm", color = wes2[3], fill = wes2[3])+
  labs(x = "Absolute Latitude (º)",
       y = "Thermal Breadth (ºT)")+
  theme_few()
breadth_lat
ggsave(filename = "breadth_lat.png", width = 30, height = 20, dpi = 300, units = "cm")

breadth_order_lat <- ggplot(data = thermal_traits_order, aes(abs(lat), thermal_breadth))+
  geom_point(aes(color = order), alpha =0.1)+
  facet_wrap(.~order)+
  geom_smooth(aes(color = order, fill = order), method = "lm")+
  theme_few()+
  scale_color_brewer(palette = "Set2")+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none")+
  labs(x = "Absolute Latitude (º)",
       y = "Thermal Breadth (ºT)")
breadth_order_lat
ggsave(filename = "breadth_order_lat.png", width = 30, height = 20, dpi = 300, units = "cm")

# 6. Plot traits variability differences  ----
# ~~ a) order ----
boxplot_traits_ord<- ggplot(data=thermal_traits_order, aes(x=order))+
  geom_boxplot(aes(y=tmin),
               fill= wes4[5],
               color = "lightgrey")+
  geom_boxplot(aes(y=tmax),
               fill= wes2[4],
               color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = "Order", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
boxplot_traits_ord
ggsave("boxplot_traits_ord.png", width = 30, height = 30, units = "cm")

# ~~ b) feeding guild ----
boxplot_traits_fg <- ggplot(data=thermal_traits_fg, aes(x=feeding_guild))+
  geom_boxplot(aes(y=tmin),
               fill= wes2[5],
               color = "lightgrey")+
  geom_boxplot(aes(y=tmax),
               fill= wes2[4],
               color = "lightgrey")+
  labs(title = "Thermal traits across taxa",
       x = "Feeding guild", 
       y = "Thermal traits (ºC)")+
  ggthemes::theme_few()+
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust=1))+
  theme(legend.position = "none")
boxplot_traits_fg
ggsave("boxplot_traits_fg.png", width = 30, height = 30, units = "cm")
