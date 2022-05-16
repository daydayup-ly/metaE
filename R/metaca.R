
#' Title  do meta-analysis by fitting random effect module
#'
#' @param data a data.frame which contain group information and values for statistic.
#' @param study a character string specifying the col.name of data, and the column contains data set information, default='Study'.
#' @param group a character string specifying the col.name of data, and the column contains group information, default='Group'
#' @param Control a character string specifying the control group.
#' @param Case   a character string specifying the case group.
#' @param value  a character string specifying the col.name of data which the column is the value for statistic,default='Value'.
#' @param ci confidence level of the interval, default=0.95.
#'
#' @return res a data.frame
#' @export
#'
#' @examples metaca(data,study='Study',group='disease_stage',Control='healthy',Case='cancer',value='diversity')
metaca<- function(data,study='Study',group='Group',Control='control',Case='case',value='Value',ci=0.95){


  # prepare data
  names(data)[names(data)==study] <- 'Study'
  names(data)[names(data)==group] <- 'Group'
  names(data)[names(data)==value] <- 'Value'

  data[which(data['Group']== Control),'Group'] <-'control'
  data[which(data['Group']== Case),'Group'] <-'case'

  data <- data %>% select(Study, Group, Value)

  # t test for each study bwt group
  data <- data %>% left_join(data %>% group_by(Study,Group) %>% summarize(mean=mean(log2(Value)))%>%
                               spread(key=Group,value=mean) %>%
                               rename(mean_log2case=case, mean_log2control=control)) %>%
    mutate(log2FC=log2(Value)-mean_log2control)

  stat <- data %>% group_by(Study) %>%
    do(broom::tidy(t.test(log2FC~Group, data=., conf.int=TRUE, conf.level=ci)))%>%
    select(Study,
           log2FC=estimate,
           Pvalue=p.value,
           CI_low=conf.low,
           CI_high=conf.high)

  # meta-analysis by fitting random effect module
  REM <- tibble(Study=character(0),
                log2FC=numeric(0), Pvalue=numeric(0),
                CI_low=numeric(0), CI_high=numeric(0))


  data$Group <- factor(data$Group,levels=(c('control','case')))
  fit<-lmerTest:::lmer(log2FC~Group+(1|Study), data)  #fit module
  cf<-confint(fit,level = ci)

  REM<-bind_rows(REM, tibble(
    Study="Meta",

    log2FC=summary(fit)$coefficients["Groupcase", "Estimate"],
    Pvalue=anova(fit)$`Pr(>F)`,
    CI_low=cf["Groupcase",1],
    CI_high=cf["Groupcase",2]
  ))

  res <- stat %>% bind_rows(REM) %>%
    mutate(Significance=case_when(
      Pvalue<0.05 & log2FC>0 ~ sprintf("* higher in %s",Case),
      Pvalue<0.05 & log2FC<0 ~ sprintf("* higher in %s",Control),
      TRUE~"ns"
    ))

  return(res)
}
