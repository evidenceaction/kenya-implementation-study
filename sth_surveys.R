library(foreach)
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lme4)
library(plm)
library(lmtest)
library(car)
library(quantreg)

source("cluster2.R")

# Prep data ---------------------------------------------------------------

sth.infection.types <- c("as", "tr", "hk", "sm", "sth")
infection.types <- c(sth.infection.types, "shaem")

prepost.sth.data <- read.csv("~/Data/Kenya STH/Y1Y2_60_prepost.csv", as.is=TRUE) %>%
  set_names(gsub("_", ".", names(.))) %>%
  rename(c("asc.high"="as.high",
           "id"="origin.id")) %>%
  mutate(id=paste(origin.id, n.survey, sep="-"),
         schoolcode=factor(schoolcode),
         districtcode=factor(districtcode),
         date=as.Date(date, format="%m/%d/%Y"),
         sex=factor(gender, levels=c("female", "male")),
         survey=factor(survey) %>% relevel(ref="Y1pre"),
         county=factor(county),
         n.survey.fac=factor(n.survey)) %>%
  mutate(shaem.high=shaeme > 50,
         shaemepg=NA)

prepost.sth.data[, paste0(sth.infection.types, "e")] <- NA

prepost.sth.data$gender <- NULL

prevalence.col <- grep("infect$", names(prepost.sth.data), value=TRUE)
prepost.sth.data[, prevalence.col] <- llply(prepost.sth.data[, prevalence.col], levels=c("positive", "negative"), factor)
prepost.sth.data[, paste0(prevalence.col, ".bin")] <- ifelse(prepost.sth.data[, prevalence.col] == "positive", 1, 0)

hi.intensity.col <- grep("high$", names(prepost.sth.data), value=TRUE)
prepost.sth.data[, hi.intensity.col] <- prepost.sth.data[, hi.intensity.col] == 1
prepost.sth.data[, paste0(hi.intensity.col, ".bin")] <- ifelse(prepost.sth.data[, hi.intensity.col] == TRUE, 1, 0)


infection.levels <- data.frame(infection.type=sth.infection.types %>% setdiff("sth"),
                               moderate.epg=c(5000, 1000, 2000, 100),
                               high.epg=c(50000, 10000, 4000, 400))


prepost.sth.long.data <- reshape(prepost.sth.data,  
                                 direction="long", 
                                 varying=list(paste0(infection.types, "infect"),
                                              paste0(infection.types, "epg"),
                                              paste0(infection.types, "e"),
                                              paste0(infection.types, ".high"),
                                              paste0(infection.types, "infect.bin"),
                                              paste0(infection.types, ".high.bin")),
                                 v.names=c("infect", "epg", "e", "high.intensity", "infect.bin", "high.intensity.bin"),
                                 times=infection.types,
                                 timevar="infection.type",
                                 ids=id) %>%
  merge(infection.levels, by="infection.type", all.x=TRUE, all.y=TRUE) %>%
  mutate(check.mod.hi=(infection.type == "sth") | ((epg >= moderate.epg) == high.intensity),
         high.intensity.only=epg >= high.epg,
         high.intensity.only.bin=ifelse(high.intensity.only == TRUE, 1, 0))

# High frequency data -----------------------------------------------------

hf.sth.data <- read.csv("~/Data/Kenya STH/Y1_HF_noschoolname.csv", as.is=TRUE) %>%
  set_names(gsub("_", ".", names(.))) %>%
  rename(c("asc.high"="as.high",
           "id"="origin.id")) %>%
  mutate(schoolcode=factor(schoolcode),
         districtcode=factor(districtcode),
         prov.code=factor(prov.code),
         dateofsurvey=as.Date(dateofsurvey, format="%m/%d/%Y"),
         treatdate=as.Date(treatdate, format="%m/%d/%Y"),
         sincetreat=as.numeric(dateofsurvey - treatdate),
         sincetreat2=sincetreat^2,
         sex=factor(gender, levels=c(2, 1), labels=c("female", "male"))) %>%
  group_by(origin.id) %>%
  mutate(num.survey=n()) %>%
  ungroup

prevalence.col <- grepl("infect$", names(hf.sth.data))
hf.sth.data[, prevalence.col] <- llply(hf.sth.data[, prevalence.col], equals, y=1)

# Plots -------------------------------------------------------------------

ggplot(prepost.sth.long.data) + 
  geom_jitter(aes(x=factor(n.survey), y=epg), alpha=0.5) +
  facet_wrap(~ infection.type) + labs(x="Survey")

prepost.sth.long.data %>% 
  filter(infection.type != "sth") %>%
  ggplot() +
  geom_jitter(aes(y=epg, x=n.survey.fac, alpha=0.15, color=factor(n.survey))) +
  geom_hline(aes(yintercept=moderate.epg), linetype="dotted", data=infection.levels) +
  labs(x="Survey", y="EPG") +
  scale_color_discrete("Survey") +
  facet_grid(infection.type ~ county, margins=TRUE) 

ggplot(prepost.sth.long.data %>% filter(epg > 10000)) + 
  geom_jitter(aes(x=factor(n.survey), y=epg, color=county == "BUSIA"), alpha=0.5) +
  facet_wrap(~ infection.type) + labs(x="Survey")

ggplot(prepost.sth.data) +
  geom_density(aes(x=asepg, color=factor(n.survey)))

ggplot(prepost.sth.data %>% filter(asepg > 10000)) +
  geom_density(aes(x=asepg, color=factor(n.survey))) +
  facet_wrap(~ county)

ggplot(prepost.sth.data %>% filter(asepg > 10000, county == "BUSIA")) +
  geom_density(aes(x=asepg, color=factor(n.survey)))

ggplot(prepost.sth.long.data) + 
  geom_bar(aes(x=factor(n.survey), fill=infect), position="fill") +
  facet_wrap(~ infection.type) + labs(x="Survey")

ggplot(prepost.sth.long.data) + 
  geom_bar(aes(x=high.intensity, fill=factor(n.survey)), position="dodge") +
#   geom_bar(aes(x=factor(n.survey), fill=high.intensity), position="fill") +
  facet_wrap(~ infection.type) + labs(x="Survey") 

ggplot(hf.sth.data %>% filter(sincetreat > 0)) +
#   geom_vline(xintercept=0) + 
  geom_jitter(aes(sincetreat, asepg)) +
  geom_smooth(aes(sincetreat, asepg)) 
#   facet_wrap(~ schoolcode)

ggplot(hf.sth.data %>% filter(sincetreat > 0)) +
#   geom_vline(xintercept=0) + 
  geom_smooth(aes(sincetreat, ifelse(asinfect, 1, 0))) 
#   facet_wrap(~ schoolcode)

# Analysis -------------------------------------------------------------

filter(prepost.sth.long.data, infection.type == "as", n.survey == 1, county == "BUSIA") %>%  
  lmer(infect ~ (1|schoolcode), data=., REML=FALSE) %>% 
  summary

for (it in infection.types) {
  sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it) 
  reg.res <- lm(infect.bin ~ n.survey.fac, data=reg.data)
  reg.res %>% coeftest(., vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
}

# epg.quant.res <- foreach (it="as") %do% { #infection.types) %do% {
#   sprintf("%s\n", it) %>% cat
  quant.reg.res <- prepost.sth.long.data %>% 
    filter(infection.type == it) %>%
#     mutate(censor=24024,
#            lcensor=0,
#            epg=pmin(censor, epg)) %>%
#     crq(Curv(epg, censor, ctype="right") ~ n.survey.fac, taus=80:99/100, data=.) 
#     crq(Curv(epg, lcensor, ctype="left") ~ survey, tau=0.92, data=., method="Pow", start="global") 
    rq(epg ~ n.survey.fac, tau=c(80, 81, 83, 86, 88, 90, 93, 95, 96, 98)/100, data=.) # %>%
#     summary(se="boot") %>%
#     print
# }

# for (it in "as") { # infection.types) {
#   sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it) 
  reg.res <- lm(epg ~ n.survey.fac, data=reg.data)
  reg.vcov <- vcov.cluster(reg.data, reg.res, cluster1="districtcode")
  reg.res %>% coeftest(., vcov=reg.vcov) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
# }

quant.predict.data <- quant.reg.res %>% 
  predict(data.frame(n.survey.fac=factor(1:4))) %>% 
  t %>% 
  as.data.frame %>% 
  mutate(tau=rownames(.), 
         tau=sub("tau=\\s+", "", tau) %>% as.numeric) %>% 
  set_names(c(paste0("quant.", 1:4), "tau")) %>% 
  reshape(direction="long", varying=1:4, idvar="tau", timevar="survey") %>% 
  filter(tau %in% c(0.9, 0.95, 0.96, 0.98))

ols.predict.data <- reg.res %>% 
  predict(data.frame(n.survey.fac=factor(1:4))) %>% 
  as.data.frame %>% 
  set_names("ave") %>% 
  mutate(survey=1:4)

ggplot(quant.predict.data, aes(x=survey, y=quant, color=factor(tau))) + 
  geom_point(aes(shape="Quant")) + 
  geom_line(aes(linetype="Quant")) +
  geom_point(aes(y=ave, shape="Mean"), data=ols.predict.data) + 
  geom_line(aes(y=ave, linetype="Mean"), data=ols.predict.data)
              
epg.quant.res <- foreach (it="as") %do% { #infection.types) %do% {
  sprintf("%s\n", it) %>% cat
  prepost.sth.long.data %>% 
    filter(infection.type == it, infect == "positive") %>%
    rq(epg ~ n.survey.fac, data=., tau=seq(0.05, 0.95, 0.05)) %>% # tau=c(50, 55, 60, 65, 75)/100) %>%
#     summary(se="boot") %>%
    summary(se="boot") %>%
    print
}

for (it in "as") { # infection.types) {
  sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it, infect == "positive") 
  reg.res <- lm(epg ~ n.survey.fac, data=reg.data)
  reg.res %>% coeftest(., vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
}

for (it in infection.types) {
  sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it) 
  reg.res <- lm(high.intensity.bin ~ n.survey.fac, data=reg.data)
  reg.res %>% coeftest(., vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
}

for (it in infection.types) {
  sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it, infect == "positive") 
  reg.res <- lm(high.intensity.bin ~ n.survey.fac, data=reg.data)
  reg.res %>% coeftest(., vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
}
