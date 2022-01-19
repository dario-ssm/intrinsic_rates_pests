#### MODELOS MIXTOS GENERALIZADOS PARA METAANÁLISIS ####

#### enfoque con nls() ####
briere_1 <- function(a, temp, Tmin, Tmax){
  a*temp*(temp-Tmin)*(Tmax-temp)^(1/2)}
# testing the equation works: 3^2 + 1 = 10



fitted_br_1 <- nls(r_ex ~ briere1(a, temp_ex, Tmin,Tmax),
                   data = example1,
                   start = list(a = 1.890e-04,
                                Tmin =5.2700,
                                Tmax= 37.670), 
                   trace = T)

summary(fitted_br_1) #problema: los tmin y tmax no ayudan... ¿por qué?





# funciones para sacar iniciales?
#initVals <- getInitial(briere_1(a,temp_ex,r_ex,Tmin,Tmax),data=example1)
#init_briere <- getInitial(r ~ briere_1(a, Temp,Tmin,Tmax),
 #                         data = data.frame(r = r_ex, Temp = temp_ex))


#Formula: rT ~ aa * T * (T - Tmin) * (Tmax - T)^(1/2)
#### devRate ####
# buscar los parámetros con el paquete devRate
library(devRate)

#buscar qué modelo es el mejor para la familia:
devRateFind(familySP =  "Noctuidae") #nos salen todas las ocurrencias
devRateInfo(eq = briere1_99) #parámetros para todos
 #vemos que para Gelechiidae, a= 2.25,Tmin=10.98,Tmax=39.93)

# para construir el modelo
nlsrbriere <- devRateModel(eq = briere1_99,
                           temp = temp_ex, 
                           devRate = example1$r_ex,
                           startValues = list(aa = 2.25, Tmin = 10.98,Tmax = 39.93))


# para sacar gráfica
devRatePlot(eq = briere1_99, nlsDR = nlsrbriere,
            temp = example1$temp_ex,
            devRate = example1$r_ex)
#para ver parámetros
devRatePrint(myNLS = nlsrbriere)
devRatePrin

# TAMBIÉN podemos usar el devRate solo para consultar initials
#vemos que para Noctuidae, a= 2.25,Tmin=10.98,Tmax=39.93)

briere_1 <- function(r, Temp,a, Tmin, Tmax){
  r = a*Temp*(Temp-Tmin)*(Tmax-Temp)^(1/2)
  df <- data.frame(Temp,r)
  print(list(r,a,Temp,Tmin,Tmax))
}

briere_1(r=r_ex,a=0.00025,Temp = temp_ex,Tmin=10.98,Tmax=39.93)

## ejemplo Ngowi et al. 2017.
rdev <- 1/c(25.55,19.30,11.78,7.89,5.51,2.71)
temp <- c(10,15,20,25,30,35)
rint <- c(0.01,0.07,0.13,0.21,0.20,0)
ngowi2017 <- data.frame(rdev,temp,rint)

ngowi_briere1_dev <- nls(rdev ~ briere1(a, temp, Tmin,Tmax),
                   data = ngowi2017,
                   start = list(a = 5e-05,
                                Tmin =8,
                                Tmax= 35), 
                   trace = T)
# sale un error de step factor. Ver esto:
library(nlme)
fit <- nls(rdev ~ briere1(a, temp, Tmin,Tmax),
           data = ngowi2017,
           start = list(a = 5e-05,
                        Tmin =8,
                        Tmax= 35))

curve(predict(fit, newdata = data.frame(Age=x)), add=TRUE)
Logistic_gnls <- gnls(rdev ~ briere_1(a, temp, Tmin,Tmax),
                      data = ngowi2017,
                      start = coef(fit))
summary(fitted_br_1)

####package devRate ####
library(devRate)
briere2NLS <- lapply(listDS, function(myDataSet){
  
  devRateModel(eq = briere1_99, 
               dfData = train_data,
               startValues = list(aa = 100, Tmin = 7.8, Tmax = 9)
  )
})
myNLS <- devRateModel(eq = briere1_99, temp = temp, devRate = growth,
                      startValues = list(aa = 0.00000002,Tmin=0.007,Tmax=0.007))



devRatePlot(eq = briere1_99, nlsDR = myNLS, temp = myT, devRate = myDev,
            spe = TRUE, pch = 16, lwd = 2, ylim = c(0, 0.10))
#### package thermPerf####
library(devtools)
#install_github("mdjbru-R-packages/thermPerf")
library(thermPerf)
models <- getModelLibrary() # This returns a list of customModel objects
myModel <- models[["briere1"]]

fits = fitModels(getModelLibrary(), temp, r)
plot(fits, xlim = c(15, 35), ylim = c(0, 1), las = 1)

f=function(temp,params){
  a=params[["a"]]
  Tmin=params[["Tmin"]]
  Tmax=params[["Tmax"]]
  return(a*temp*(temp-Tmin)*(Tmax-temp)^(1/2))
}
myModel=buildModel(f,"briere-1",
                   r~a*temp*(temp-Tmin)*(Tmax-temp)^(1/2),
                   c("a","Tmin","Tmax"),list(a=0.0002,Tmin=7.5,Tmax=34))
fitModel(customModel = myModel,x = temp,y = growth,
         initParams = NULL)