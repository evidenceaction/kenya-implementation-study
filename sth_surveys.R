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
library(yaml)

source("cluster2.R")

# Init --------------------------------------------------------------------

config <- yaml.load_file("local_config.yaml")

western.province <- c("KAKAMEGA", "VIHIGA", "BUNGOMA", "BUSIA")

# Prep data ---------------------------------------------------------------

sth.infection.types <- c("as", "tr", "hk", "sm", "sth")
infection.types <- c(sth.infection.types, "shaem")

prepost.sth.data <- read.csv(sprintf("%s/Kenya STH/Y1Y2_60_prepost.csv", config$data_path), as.is=TRUE) %>%
  set_names(gsub("_", ".", names(.))) %>%
  rename(as.high=asc.high,
         origin.id=id) %>%
  mutate(id=paste(origin.id, n.survey, sep="-"),
         schoolcode=factor(schoolcode),
         districtcode=factor(districtcode),
         date=as.Date(date, format="%m/%d/%Y"),
         survey.mon=format(date, "%Y-%m") %>% factor(),
         sex=factor(gender, levels=c("female", "male")),
         survey=factor(survey, levels=c("Y1pre", "Y1post", "Y2pre", "Y2post")), # relevel(ref="Y1pre"),
         county=factor(county),
         n.survey.fac=factor(n.survey),
         deworm.status=sub("Y\\d", "", as.character(survey)) %>% factor %>% relevel(ref="pre")) %>%
  mutate(shaem.high=shaeme > 50,
         shaemepg=NA) %>%
  select(-matches(sprintf("(%s)e$", paste(infection.types, collapse="|")))) %>%
  select(-gender) 

prevalence.col <- grep("infect$", names(prepost.sth.data), value=TRUE)
prepost.sth.data[, prevalence.col] <- llply(prepost.sth.data[, prevalence.col], levels=c("positive", "negative"), factor)
prepost.sth.data[, paste0(prevalence.col, ".bin")] <- ifelse(prepost.sth.data[, prevalence.col] == "positive", 1, 0)

hi.intensity.col <- grep("high$", names(prepost.sth.data), value=TRUE)
prepost.sth.data[, hi.intensity.col] <- prepost.sth.data[, hi.intensity.col] == 1
prepost.sth.data[, paste0(hi.intensity.col, ".bin")] <- ifelse(prepost.sth.data[, hi.intensity.col] == TRUE, 1, 0)

infection.levels <- data.frame(infection.type=sth.infection.types %>% setdiff("sth"),
                               moderate.epg=c(5000, 1000, 2000, 100),
                               high.epg=c(50000, 10000, 4000, 400))

prepost.sth.long.data <- prepost.sth.data %>% 
  (l(.data ~ reshape(.data, 
                       direction="long",
                       varying=llply(.v.names, function(vn, it) paste0(it, vn), infection.types),
                       v.names=.v.names %>% gsub("^\\.", "", .), 
                       times=infection.types,
                       timevar="infection.type",
                       ids=.data$id)) %>% add_args(.v.names=c("infect", "epg", ".high", "infect.bin", ".high.bin"))) %>%
  merge(infection.levels, by="infection.type", all.x=TRUE) 

prepost.school.data <- prepost.sth.long.data %>%
  group_by(schoolcode, n.survey, infection.type) %>%
  summarize(mean.epg=mean(epg), mean.infect.bin=mean(infect.bin))

# High frequency data -----------------------------------------------------

hf.sth.data <- read.csv(sprintf("%s/Kenya STH/Y1_HF_noschoolname.csv", config$data_path), as.is=TRUE) %>%
  set_names(gsub("_", ".", names(.))) %>%
  rename(as.high=asc.high,
         origin.id=id) %>%
  mutate(schoolcode=factor(schoolcode),
         districtcode=factor(districtcode),
         prov.code=factor(prov.code),
         dateofsurvey=as.Date(dateofsurvey, format="%m/%d/%Y"),
         treatdate=as.Date(treatdate, format="%m/%d/%Y"),
         sincetreat=as.numeric(dateofsurvey - treatdate),
         sincetreat2=sincetreat^2,
         sex=factor(gender, levels=c(2, 1), labels=c("female", "male")),
         deworm.status=factor(sincetreat >= 0, levels=c(TRUE, FALSE), labels=c("post", "pre"))) %>%
  select(-gender) %>%
  group_by(origin.id) %>%
  mutate(num.survey=n()) %>%
  ungroup

prevalence.col <- grepl("infect$", names(hf.sth.data))
hf.sth.data[, prevalence.col] <- llply(hf.sth.data[, prevalence.col], equals, y=1)

hf.sth.long.data <- hf.sth.data %>%
  mutate(shaemepg=NA, shaem.high=NA, sth.high=NA) %>%
  (l(.data ~ reshape(.data,
                     direction="long",
                     varying=llply(.v.names, function(vn, it) paste0(it, vn), infection.types),
                     v.names=.v.names %>% sub("^\\.", "", .),
                     ids=c("origin.id", "dateofsurvey"),
                     times=infection.types,
                     timevar="infection.type",
                     new.row.names=rep(rownames(.data), each=length(infection.types)) %>% paste(infection.types, sep="."))) %>% 
     add_args(.v.names=c("infect", "epg", ".high"))) %>%
  merge(infection.levels, by="infection.type", all.x=TRUE) %>%
  mutate(high=epg > moderate.epg) %>%
  merge(filter(., deworm.status == "pre") %>% 
          arrange(nr.survey) %>% 
          group_by(origin.id, infection.type) %>% 
          summarize(bl.epg=last(epg), 
                    bl.high=last(high)),
        by=c("origin.id", "infection.type"),
        all.x=TRUE)

# Save data ---------------------------------------------------------------

save(list=ls(pattern="^(prepost|hf).+data$"), file="cleaned_sth.RData")

# Plots -------------------------------------------------------------------

prepost.sth.long.data %>%
  ggplot()

ggplot(prepost.sth.long.data) + 
  geom_jitter(aes(x=factor(n.survey), y=epg), alpha=0.5) +
  facet_wrap(~ infection.type) + labs(x="Survey")

prepost.sth.long.data %>% 
  filter(!infection.type %in% c("sth", "shaem"),
         county %in% western.province) %>%
  ggplot() +
  geom_jitter(aes(y=epg, x=n.survey.fac, color=factor(deworm.status)), alpha=0.4) +
  geom_hline(aes(yintercept=moderate.epg), linetype="dotted", data=infection.levels) +
  labs(x="Survey", y="EPG") +
  scale_color_discrete("Survey") +
  facet_grid(~ infection.type, margins=TRUE) 

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

filter(prepost.sth.long.data, infection.type == "sth", n.survey == 1) %>%  
  filter(county %in% western.province) %>%
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
    filter(infection.type == "as") %>%
#     mutate(censor=24024,
#            lcensor=0,
#            epg=pmin(censor, epg)) %>%
#     crq(Curv(epg, censor, ctype="right") ~ n.survey.fac, taus=80:99/100, data=.) 
#     crq(Curv(epg, lcensor, ctype="left") ~ survey, tau=0.92, data=., method="Pow", start="global") 
    rq(epg ~ n.survey.fac, tau=c(80, 81, 83, 86, 88, 90, 93, 95, 96, 98)/100, data=.) # %>%
#     summary(se="boot") %>%
#     print
# }

for (it in "as") { # infection.types) {
#   sprintf("%s\n", it) %>% cat
  reg.data <- prepost.sth.long.data %>% filter(infection.type == it) 
  reg.res <- lm(epg ~ n.survey.fac, data=reg.data)
  reg.vcov <- vcov.cluster(reg.data, reg.res, cluster1="districtcode")
  reg.res %>% coeftest(., vcov=reg.vcov) %>% print
  reg.res %>% lht("n.survey.fac3 = n.survey.fac4", vcov=vcov.cluster(reg.data, ., cluster1="districtcode")) %>% print
}

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

for (it in "sth") { # infection.types) {
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

reg.data <- hf.sth.long.data %>% 
  filter(infection.type == "as", sincetreat >= 0) %>% 
  mutate(sincetreat_bl.epg=sincetreat * bl.epg,
         sincetreat_bl.high=sincetreat * bl.high)

predict.data <- expand.grid(dateofsurvey=mean(reg.data$dateofsurvey), 
                            sincetreat=range(reg.data$sincetreat),
                            bl.high=0:1) %>%
#                             bl.epg=rep(quantile(reg.data$bl.epg, seq(0, 1, 0.05)), each=2)) %>% 
            mutate(#sincetreat_bl.epg=sincetreat * bl.epg,
                   sincetreat_bl.high=sincetreat * bl.high) 
predicted.epg <- lm(epg ~ dateofsurvey + sincetreat + sincetreat_bl.high, reg.data) %T>% 
  (l(reg.res ~ coeftest(reg.res, vcov=vcov.cluster(reg.data, reg.res, cluster1="districtcode")) %>% print)) %T>%
  (l(reg.res ~ lht(reg.res, c("sincetreat_bl.high", "sincetreat")) %>% print)) %>%
  predict(newdata=predict.data) %>% unname

predict.data %>%
  mutate(epg=predicted.epg) %>% 
  ggplot(aes(x=sincetreat, y=epg)) + geom_line(aes(color=factor(bl.high))) 

qplot(x=predict.data$sincetreat, y=., geom="point")
  ggplot(predict.data, y=.) %>% geom_line(aes(x=sincetreat, color=bl.epg))

d_ply(prepost.sth.long.data, .(infection.type), function(df) {
  lm(infect.bin ~ deworm.status + survey.mon, data=df) %>% 
    coeftest(vcov=vcov.cluster(df, ., cluster1="districtcode")) %>%
    print
})  

prepost.sth.long.data %>% 
  filter(county %in% western.province, 
         infection.type %in% c("sth", "as")) %>%
  d_ply(.(infection.type), function(df) {
    cat(sprintf("%s\n", df$infection.type[1]))
    lm(high.bin ~ deworm.status + survey.mon, data=df) %>% 
      (l(reg.res ~ {
        robust.vcov <- vcov.cluster(df, reg.res, cluster1="districtcode")
        coeftest(reg.res, vcov=robust.vcov) %>% print 
#         reg.res %>% lht(grep("survey\\.mon", coef(.) %>% names, value=TRUE), vcov=robust.vcov) %>% print
      })) 
  })

d_ply(prepost.sth.long.data, .(infection.type), function(df) {
  try(rq(epg ~ deworm.status + survey.mon, data=df, tau=c(0.80, 0.98)) %>% summary %>% print)
})  

lm(infect.bin ~ deworm.status + survey.mon, data=prepost.sth.long.data, subset=infection.type == "as")

# hf.sth.long.data %>%
#   filter(sincetreat %in% c(90:120, 270:300)) %>%
#   mutate(treat=sincetreat <= 120,
#          high.bin=ifelse(high == TRUE, 1, 0)) %>%
#   lm(high.bin ~ treat, data=.)