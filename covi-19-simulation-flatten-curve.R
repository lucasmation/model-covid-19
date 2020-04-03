#Lucas Mation, 2020-mar

install.packages("EpiModel", dependencies = TRUE)
library(EpiModel)
library(tidyverse)
library(data.table)


#documentation
#This similatior is just copied from this blog post
#Great blog post: https://timchurches.github.io/blog/posts/2020-03-10-modelling-the-effects-of-public-health-interventions-on-covid-19-transmission-part-1/
#another great one


#Ventilator need per infected calculation

#1)https://sccm.org/Blog/March-2020/United-States-Resource-Availability-for-COVID-19
# A recent AHA webinar on COVID-19 projected that 30% (96 million) of the U.S. population will test positive, 
# with 5% (4.8 million) being hospitalized. Of the hospitalized patients, 40% (1.9 million) would be admitted 
# to the ICU, and 50% of the ICU admissions (960,000) would require ventilatory support.
# 0.3*0.05*0.4*0.5 = 0.003

#2) WHO-china mission, 
# https://www.who.int/publications-detail/report-of-the-who-china-joint-mission-on-coronavirus-disease-2019-(covid-19)
#Most people infected with COVID-19 virus have mild disease and recover. Approximately 80% of laboratory confirmed 
# patients have had mild to moderate disease, which includes non-pneumonia and pneumonia cases, 13.8% have severe 
# disease (dyspnea, respiratory frequency ≥30/minute, blood oxygen saturation ≤93%, PaO2/FiO2 ratio <300, and/or
#lung infiltrates >50% of the lung field within 24-48 hours) and 6.1% are critical (respiratory
                                                                                                                                            failure, septic shock, and/or multiple organ dysfunction/failure)
#3) Imperial College reports
#https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/
#https://www.imperial.ac.uk/news/196496/coronavirus-pandemic-could-have-caused-40/
#

#Flatten the curve plotsÇ
#https://static01.nyt.com/images/2020/03/22/science/11SCI-VIRUS-CURVE1/11SCI-VIRUS-TRACKER1-jumbo.jpg?quality=90&auto=webp




control <- control.icm(type = "SIR"  ,  #model
                       nsteps = 200,    #n periods
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
                   di.rate = 0.01/13,#deaths among infected
                   dr.rate = 0    #deaths among recorvered
                   )  

sim <- icm(param, init, control)
plot(sim)

sim_df <- sim[[3]] %>% bind_cols() %>% data.table
names(sim_df) <- c('s.num','i.num','num','r.num','si.flow','ir.flow','ds.flow','di.flow','dr.flow','a.flow')
sim_df[,sum(di.flow)]
sim_df[,num.br:=num*210]
sim_df[,i.num2.br:=i.num*210]
sim_df[,i.hos.br:=i.num2.br*0.146]
sim_df[,i.icu.br:=i.num2.br*0.025]
sim_df[,t:=1:.N]


sim_df %>% ggplot(aes(t,i.icu.br)) + geom_line(color='red', size=1) + 
  geom_hline(yintercept = c(30),size=1) + 
    theme_light() + ylab('milhares') + xlab('dias') +
    annotate("text", x = 70, y = 2500, label = "Demanda por UTI sem medidas de",col="red",size=5.5) + 
    annotate("text", x = 70, y = 2300, label = "contenção (2.5% dos infectados)",col="red",size=5.5) + 
    annotate("text", x = 50, y = 160, label = "Capacidade do sistema",col="black",size=5.5) +
  labs(caption='Brasil, Modelo SIR (R pacote EpiModel), R_0=2.24, supondo 30mil leitos de UTI',size=5) +
  
#Obs, em log os flatten the curve plots ate que fazem sentido.   scale_y_continuous(trans='log')