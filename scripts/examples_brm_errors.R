library(nlme)
library(brms) # at least 0.8
library(ggplot2)

f1 <- circumference ~ phi1 / (1 + exp(-(age - phi2)/phi3))

prior_1 <- c(set_prior("normal(200, 50)", nlpar = "phi1"),
             set_prior("normal(700, 50)", nlpar = "phi2"),
             set_prior("normal(350, 50)", nlpar = "phi3"))
n1_b <- brm(f1, 
            data = Orange,
            nonlinear = phi1 + phi2 + phi3 ~ 1,
            prior = prior_1, 
            chains = 3)


2 <- nlme(f1,
          data = Orange,
          fixed = phi1 + phi2 + phi3 ~ 1,
          random = phi1 ~ 1,
          groups = ~ Tree,
          start = coef(n1))

summary(n2)

# prior_1 plus prior for random effect, is there a way to update a  `prior_frame` ?

prior_2 <- c(set_prior("normal(200, 50)", nlpar = "phi1"),
             set_prior("normal(700, 50)", nlpar = "phi2"),
             set_prior("normal(350, 50)", nlpar = "phi3"),
             set_prior("cauchy(30,2)", class = "sd"))
n2_b <- brm(bf(f1, nonlinear = list(phi1 ~ (1|Tree),
                                    phi2 ~ 1,
                                    phi3 ~ 1),nl=TRUE),
            data = Orange,
            prior = prior_2,
            chains = 3,
            iter = 4000)
summary(n2_b)


#### check other option with Orange Example ####

f1 <- circumference ~ phi1 / (1 + exp(-(age - phi2)/phi3))

prior_1 <- c(set_prior("normal(200, 50)", nlpar = "phi1"),
             set_prior("normal(700, 50)", nlpar = "phi2"),
             set_prior("normal(350, 50)", nlpar = "phi3"))

orange_formula <- bf(circumference ~ phi1 / (1 + exp(-(age - phi2)/phi3)),
                             # Nonlinear variables
                              phi1 + phi2 + phi3 ~ 1,
                             # Nonlinear fit
                              nl = TRUE)

bayes_fit <- brm(
  orange_formula,
  family=gaussian(), 
  data = Orange,
  prior = prior_1 )
summary(bayes_fit) # works :____)

Rejecting initial value:
  Chain 1:   Error evaluating the log probability at the initial value.
Chain 1: Exception: normal_lpdf: Location parameter[1] is nan, but must be finite!  (in 'model1e6847002cd7_d154ae6bed8a90c4055e64601dfe70e8' at line 42)





abracadabra <- raster::raster("/Users/Ecologia/Desktop/DARÍO_actualizada septiembre 2021/phenoInsects_brassica/spatialmaps/past/gens_total_present_geotiff.asc")
abracadabra
map_present_voltinism <- gplot(abracadabra) + 
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = c("white",rev(heat.colors(9))))+
  coord_equal()+
  theme_minimal()+
  labs(x="lon", y="lat", title = "Voltinism (Present)",
       subtitle= "Plutella xylostella, Spain",fill="Number of generations")
map_present_voltinism

writeRaster(abracadabra, "alakazam",
            format = "GTiff",
            overwrite=TRUE)
