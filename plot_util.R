library(plyr)
library(dplyr)

plot.sth.incidence <- function(.data, target.counties=NULL) {
  .data %>% 
    filter(infection.type %in% c("sth", "as", "hk")) %>%
    (l(.data ~ {
      if (!is.null(target.counties)) {
        filter(.data, county %in% western.province)
      } else {
        .data
      }
    })) %>% 
    ddply(.(infection.type), function(df) {
      predict.data <- expand.grid(survey=levels(df$survey)) 
      
      rbind(lm(high.bin ~ survey, data=df) %>% 
          (l(reg.res ~ {
            robust.vcov <- vcov.cluster(df, reg.res, cluster1="districtcode")
            
    #         lht(reg.res, "survey[T.Y2pre] = survey[T.Y2post]", vcov=robust.vcov) %>% print
            
            predict.rob(reg.res, vcov=robust.vcov, newdata=predict.data) %>% as.data.frame %>% cbind(predict.data) %>% mutate(dep.var="high.intensity")
          })), 
        lm(infect.bin ~ survey, data=df) %>% 
          (l(reg.res ~ {
            robust.vcov <- vcov.cluster(df, reg.res, cluster1="districtcode")
            
    #         lht(reg.res, "survey[T.Y2pre] = survey[T.Y2post]", vcov=robust.vcov) %>% print
            
            predict.rob(reg.res, vcov=robust.vcov, newdata=predict.data) %>% as.data.frame %>% cbind(predict.data) %>% mutate(dep.var="incidence")
          })))
    }) %>%
    mutate(infection.type=factor(infection.type, levels=c("as", "hk", "sth"), labels=c("Ascaris", "Hookworm", "All"))) %>%
    (l(ready.data ~ {
      dodge <- position_dodge(width=0.5)
      ggplot(ready.data, aes(x=survey, group=dep.var)) + 
        geom_bar(aes(y=fit, fill=dep.var), stat="identity", position=dodge, width=0.5) +
        geom_errorbar(aes(ymin=fit.min, ymax=fit.max), position=dodge, width=0.3) +
        scale_y_continuous("Proportion") + 
        scale_x_discrete("Survey") + 
        scale_fill_discrete("Proportion of", labels=c("Moderate or High Intensity Infection", "Any Infection (Prevalence)")) +
  #       facet_grid(dep.var ~ infection.type) 
        facet_wrap(~ infection.type) 
    }))
}

plot.sth.quant <- function(.data, target.counties=NULL, plot.infection.types=NULL, moderate.intensity=NULL) {
  .data %>% 
    (l(.data ~ if (!is.null(target.counties)) filter(.data, county %in% target.counties) else .data)) %>%
    (l(.data ~ if (!is.null(plot.infection.types)) filter(.data, infection.type %in% plot.infection.types) else .data)) %>%
    (l(.data ~ {
      rq(epg ~ survey, tau=tau.range, data=.data) %>%
      predict(data.frame(survey=levels(prepost.sth.long.data$survey))) %>% 
      t %>% 
      as.data.frame %>% 
      mutate(tau=rownames(.), 
             tau=sub("tau=\\s+", "", tau) %>% as.numeric) %>% 
      set_names(c(paste0("quant.", 1:4), "tau")) %>% 
      reshape(direction="long", varying=1:4, idvar="tau", timevar="survey") %>% 
      filter(tau %in% tau.range) %>%
      ggplot(aes(x=survey, y=quant, color=factor(tau))) + 
      geom_point() + 
      geom_line() +
      scale_x_discrete("Survey", labels=levels(prepost.sth.long.data$survey)) +
      scale_y_continuous("EPG") +
      scale_color_discrete("Percentile", labels=sprintf("%d%%", tau.range * 100)) 
    #   geom_point(aes(y=ave, shape="Mean"), data=ols.predict.data) + 
    #   geom_line(aes(y=ave, linetype="Mean"), data=ols.predict.data)
    }) %>% add_args(tau.range=c(80, 90, 92, 95, 98)/100)) %>%
    (l(plot.obj ~ {
      if (!is.null(moderate.intensity)) {
        plot.obj + 
          geom_hline(yintercept=moderate.intensity, linetype="dotted") +
          annotate("text", x=4, y=moderate.intensity, label="Moderate Intensity", size=4, vjust=-0.8, hjust=0.2) 
      } else {
        plot.obj
      }
    }))
}