library(foreign)
library(magrittr)
library(plyr)
library(dplyr)
library(ICC)
library(lme4)

miguel.kremer.residual <- 0.273 / sqrt(1 - 0.23)

calc.mdes <- function(sig.level=0.05, power=0.8, alloc.frac=0.5, num.clust=300, clust.size=25, icc=0.04) {
  M <- qt(sig.level/2, df=num.clust - 2, lower.tail=FALSE) + qt(power, df=num.clust - 2)

  (M / sqrt(alloc.frac * (1 - alloc.frac) * num.clust)) * sqrt(icc + (((1 - icc) / clust.size)))
}  

calc.mde <- function(sig.level=0.05, power=0.8, alloc.frac=0.5, num.clust=300, clust.size=25, icc=0.04, residual=miguel.kremer.residual) {
  M <- qnorm(sig.level/2, lower.tail=FALSE) + qnorm(power)

  (M / sqrt(alloc.frac * (1 - alloc.frac) * num.clust)) * sqrt(icc + (((1 - icc) / clust.size) * residual))
}  
 
wide.hh.data <- read.dta("~/Data/KDHS/2008/KEHR52FL.DTA", convert.underscore=TRUE, convert.factors=FALSE) %>%
  select(hhid, hv001, hv002, hv022, starts_with("hv105"), starts_with("hv121")) %>%
  rename(cluster.num=hv001, 
         hh.num=hv002,
         region=hv022)

long.hh.data <- reshape(wide.hh.data, 
                        varying=list(grep("hv105", names(wide.hh.data)),
                                     grep("hv121", names(wide.hh.data))),
                        v.names=c("age", "attended.school"),
                        direction="long") 

long.hh.data <- long.hh.data %>%
  filter(!is.na(attended.school), 
         age %in% 5:15, 
         attended.school != 9) %>% 
  arrange(id) %>%
  mutate(attended.school=ifelse(attended.school == 2, 1, 0),
         region.2=ifelse(region %in% 3:4, 1, ifelse(region %in% 13:14, 2, 3)))

# long.hh.data$attended.school <- ifelse(long.hh.data$attended.school == 2, 1, 0)

# power.params <- ICCest(cluster.num, attended.school, long.hh.data)

filter(long.hh.data, region %in% 13:14) %>%  # focusing on the Western region
  lmer(attended.school ~ (1|cluster.num), data=., REML=FALSE) %>% 
  summary

within.sd <- 0.03215
between.sd <- 0.25835

overall.sd <- filter(long.hh.data, region %in% 13:14) %>% use_series(attended.school) %>% sd

icc <- within.sd^2/(within.sd^2 + between.sd^2)

mde <- calc.mde(alloc.frac=2/3, num.clust=225, clust.size=225*0.2, 
                icc=icc, residual=between.sd)
#     icc=power.params$ICC, residual=sqrt(power.params$vara))
(mde)
(mde/overall.sd)

mdes <- calc.mdes(alloc.frac=2/3, num.clust=225, clust.size=225*0.2, icc=icc)

# filter(long.hh.data, region.2 %in% 1:2) %>%  
#   lmer(attended.school ~ (1|region.2) + (1|cluster.num), data=., REML=FALSE) %>% 
#   summary
# 
# within.sd <- 0.05941
# between.sd <- 0.2633
# 
# mde <- calc.mde(alloc.frac=2/3, num.clust=450, clust.size=225*0.2, 
#                 icc=within.sd^2/(within.sd^2 + between.sd^2), residual=between.sd)

# long.hh.data %>% 
#   lmer(attended.school ~ (1|cluster.num), data=., REML=FALSE) %>% 
#   summary
