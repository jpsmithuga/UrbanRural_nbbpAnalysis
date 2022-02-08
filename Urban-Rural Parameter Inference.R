


#Bring in likelihood functions
source("Likelihood Functions.R")

rural <- read.csv("rural_overlapping.csv")
urban <- read.csv("urban_overlapping.csv")
total <- rbind(urban, rural)


## Parameter inference
R1 <- 0.01; R2 <- 2.0; Rrange <- seq(R1, R2, by = 0.01)
k1 <- 0.01; k2 <- 1.1; krange <- seq(k1, k2, by = 0.01)

#Joint surface likelihoods
ll_total <- surflike(total, Rrange = Rrange, krange = krange)
ll_urban <- surflike(urban, Rrange = Rrange, krange = krange)
ll_rural <- surflike(rural, Rrange = Rrange, krange = krange)

#Maximum likelihood
llmax_total <- ll_total == max(ll_total)
llmax_urban <- ll_urban == max(ll_urban)
llmax_rural <- ll_rural == max(ll_rural)

# Inferred parameters
total.ests <- calc_profile(ll_total, llmax_total, conf.interval = 95, Rrange=Rrange, krange = krange)
urban.ests <- calc_profile(ll_urban, llmax_urban, conf.interval = 95, Rrange = Rrange, krange = krange)
rural.ests <- calc_profile(ll_rural, llmax_rural, conf.interval = 95, Rrange = Rrange, krange = krange)
