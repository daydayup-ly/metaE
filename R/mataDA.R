
#' Title
#'
#' @param data OTU/gene/pathway abundance table
#' @param md metadata of each samples
#' @param n feature appears in more than n data sets
#'
#' @return DA a list of data.frame of results of each data sets and meta-analysis
#' @export
#'
#' @examples
metaDA <- function(data,md,n=3){

  DA <- list()
  DA$per_study_res <- data.frame()
  study <- unique(md$Study)
  for (st in study){
    md_s <- subset(md,Study==st)
    sample_ids <- intersect(row.names(md_s),names(data))
    md_s <- md_s[sample_ids,]
    data_s <- data[,sample_ids]
    data_s <-  data_s+1

    Group <-factor(md_s$Group,levels = c('control','case'))
    group <- md_s['Group']

    format <- DESeqDataSetFromMatrix(data_s, group, design = ~Group)
    analyze <- DESeq(format)

    diff_res <- results(analyze)

    diff_res <- as.data.frame(diff_res) %>%
      mutate('diff_otu'=row.names(diff_res),'Study'=st)

    DA$per_study_res <- rbind(DA$per_study_res,diff_res)
  }


  DA$per_study_res <- DA$per_study_res %>% mutate(CI_low= DA$per_study_res$log2FoldChange-1.95*DA$per_study_res$lfcSE,
                                                  CI_hig= DA$per_study_res$log2FoldChange+1.95*DA$per_study_res$lfcSE)

  # select otu appears in more than n data sets
  r <- DA$per_study_res %>% group_by(diff_otu) %>% summarize(count=n()) %>% subset(count>=n)

  # do Meta-analysis for otu
  DA$Per_otu <- data.frame()
  otu <- c()
  for (i in r$diff_otu){
    per_otu <- DA$per_study_res %>% subset(diff_otu==i)
    ma_model_1 <- rma(log2FoldChange,lfcSE, data = per_otu)
    p <- ma_model_1$pval

    if (p < 0.01){
      otu <- c(otu,i)
      per_otu_combind <- data.frame(baseMean='', log2FoldChange=ma_model_1$b[1],
                                    lfcSE='', stat='',pvalue='',
                                    padj=p, diff_otu=i,Study ='Combined', CI_low=ma_model_1$ci.lb,CI_hig=ma_model_1$ci.ub)
      per_otu_res <- rbind(per_otu,per_otu_combind)
      DA$Per_otu <- rbind(DA$Per_otu,per_otu_res)

    }
  }

  DA$Per_otu <- DA$Per_otu %>% mutate(Significant=case_when(
    padj<=0.01 & log2FoldChange>0 ~ "case enriched",
    padj<=0.01 & log2FoldChange<0 ~ "control enriched",
    TRUE~"ns"
  ))

  return(DA)
}
