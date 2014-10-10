library(sp)
library(plyr)
library(dplyr)
library(rgeos)
library(foreach)

load("~/Data/Kenya ADM/KEN_adm1.RData")
province.adm <- gadm
rm(gadm)

load("~/Data/Kenya ADM/KEN_adm2.RData")
county.adm <- gadm
rm(gadm)

busia.id <- 46

load("~/Data/Kenya ADM/KEN_adm3.RData")
subcounty.adm <- gadm
rm(gadm)

busia.subcounties.adm <- subcounty.adm[subcounty.adm$ID_2 == busia.id, ]

load("~/Data/Kenya ADM/KEN_adm4.RData")
location.adm <- gadm 
rm(gadm) 

busia.locations.adm <- location.adm[location.adm$ID_2 == busia.id, ]

load("~/Data/Kenya ADM/KEN_adm5.RData")
sublocation.adm <- gadm
rm(gadm)

busia.sublocations.adm <- sublocation.adm[sublocation.adm$ID_2 == busia.id, ]