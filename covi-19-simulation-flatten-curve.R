######################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list with one component being a vector of length two containing
#  dS(t)/dt and dI(t)/dt
######################################################################

sir <- function(t, y, parms) {
  beta <- parms[1]
  gamma <- parms[2]
  S <- y[1]
  I <- y[2]
  return(list(c(S = -beta * S * I, I = beta * S * I - gamma * I)))
}

# Population size 
N <- 1e6 
# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma <- 1/5 
# Infectious contact rate - beta = R0/N*gamma and when R0 \approx 2.25 then  2.25/N*gamma
beta <- 4.5e-07 
# R0 for the beta and gamma values
R0 <- beta*N/gamma



suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(tidyverse))

# Grid where to evaluate
max_time <- 150
times <- seq(0, max_time, by=0.01)

# Solve ODE system using Runge-Kutta numerical method.
ode_solution <- rk4(y = c(N - 10, 10), times = times, func = sir, parms = c(beta, gamma)) %>%
  as.data.frame() %>%
  setNames(c("t", "S", "I")) %>%
  mutate(beta = beta, gama = gamma, R0 = N * beta / gamma, s = S / N, i = I / N, type = "without_intervention")






Segundo o Ministério da Saúde, país tem 14,8 mil leitos de UTIs adulto. 
Pasta abriu licitação para contratar 2.000 novos leitos, abaixo dos 2.960 que associação diz serem necessários.

stimates so far show that about 6% of people who have COVID-19 get critically sick. And about 1 in 4 of them may need a ventilator to help them breathe. But the picture is changing quickly as the infection continues to spread around the globe.





library(tidyverse)





install.packages("EpiModel", dependencies = TRUE)
library(EpiModel)

#documentation
#Great blog post: https://timchurches.github.io/blog/posts/2020-03-10-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/
#another great one


control <- control.icm(type = "SIR"  ,  #model
                       nsteps = 100,    #n periods
                       nsims = 1)      #n simulations
N <- 1000
init <- init.icm(s.num = N-1, # susceptivble 
                 i.num = 1,   # infected
                 r.num = 0)   # recorvered

param <- param.icm(inf.prob = 0.05, 
                   act.rate = 10, 
                   rec.rate = 1/20, 
                   a.rate  = 0   ,#arrival rate, no births or imigrants
                   ds.rate = 0   ,#deaths among succeptible (not infected)
                   di.rate = 0.01,#deaths among infected
                   dr.rate = 0    #deaths among recorvered
                   )  

sim <- icm(param, init, control)

icmN <- function(N){
  init <- init.icm(s.num = N-1, # susceptivble 
                   i.num = 1,   # infected
                   r.num = 0)   # recorvered
  
  param <- param.icm(inf.prob = 0.05, 
                     act.rate = 10, 
                     rec.rate = 1/20, 
                     a.rate  = 0   ,#arrival rate, no births or imigrants
                     ds.rate = 0   ,#deaths among succeptible (not infected)
                     di.rate = 0.01,#deaths among infected
                     dr.rate = 0    #deaths among recorvered
                     )  
  
  icm(param, init, control) %>% return
  
}

sim <- icmN(10000)

sim[[3]][[1]] %>% str

plot(sim)


remotes::install_github("olafmersmann/microbenchmark")
install.packages('microbenchmark')
library(microbenchmark)

microbenchmark("10^3"  =  { sim10_3 <- icmN(10^3)},
               "10^4"  =  { sim10_4 <- icmN(10^4)},
               "10^5"  =  { sim10_5 <- icmN(10^5)},
               "10^6"  =  { sim10_6 <- icmN(10^6)},
               "10^7"  =  { sim10_7 <- icmN(10^7)},
               times = 1
               )


microbenchmark("10^8"  =  { sim10_8 <- icmN(10^8)},
               times = 1
               )


plot(sim10_3)
plot(sim10_4)
plot(sim10_5)
plot(sim10_6)
plot(sim10_7)

b <- icmN(10000),
               "pseudoinverse" = {
                 b <- solve(t(X) %*% X) %*% t(X) %*% y
               },
               "linear system" = {
                 b <- solve(t(X) %*% X, t(X) %*% y)
               },
               check = check_for_equal_coefs)




