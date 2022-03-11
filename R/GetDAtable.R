#' Genetic effects analysis of miRNA: |d/a| (method 1)
#'
#' The additive (a) and dominant (d) values were calculated by the expression level of each miRNA. Edwards et al. proposed that the "|d/a|" can be used as the criterion to estimate the expression patterns of miRNAs. Specific classification criteria are as follows, |d/a| <= 0.2, additivity; |d/a| > 0.2 and |d/a| <= 0.8, partial dominance; |d/a| > 0.8 and |d/a| <= 1.2, dominance; |d/a| > 1.2, overdominance.
#'
#' @param P1_RPM A dataframe. The rpm data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the rpm of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the P2 species.
#' @param F1_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the F1 species.
#' @param rpm_threshold A numeric. the average of rpm value among all the biological replicates. By default, the average of rpm more than or equal to 1 is retained.
#'
#' @return A dataframe. The output results contain the value of "|d/a|" and grouping results for each miRNA expressed in all species (average_rpm >= rpm_threshold).
#' @export
#' @examples
#' ##Get the classification results based on the value of |d/a|
#' DAresult <- GetDAtable(P1_RPM = P1_miRNA_rpm,
#'                        P2_RPM = P2_miRNA_rpm,
#'                        F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
GetDAtable <- function(P1_RPM,P2_RPM,F1_RPM,rpm_threshold = 1){
  func_mirnafiliter <- function(data1,threshold,group){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- cbind(data2[,1],as.data.frame(apply(data2[,c(-1)],1,mean)))
    colnames(result) <- c("sequence",group)
    return(result)
  }
  P1_mirna <- func_mirnafiliter(data1 = P1_RPM,threshold = rpm_threshold,group = "P1_average_rpm")
  P2_mirna <- func_mirnafiliter(data1 = P2_RPM,threshold = rpm_threshold,group = "P2_average_rpm")
  F1_mirna <- func_mirnafiliter(data1 = F1_RPM,threshold = rpm_threshold,group = "F1_average_rpm")
  ####################################
  Parent <- merge.data.frame(P1_mirna,P2_mirna,all=TRUE)
  filter_result <- merge.data.frame(Parent,F1_mirna,all=TRUE)
  fina_result <- na.omit(filter_result)
  #####
  fina_result$d <- fina_result$F1_average-((fina_result$P1_average+fina_result$P2_average)/2)
  fina_result$a <- (fina_result$P1_average-fina_result$P2_average)/2
  fina_result$`abs(d/a)` <- abs(fina_result$d/fina_result$a)
  da <- fina_result$`abs(d/a)`
  if (length(da[da<=0.2]>0)){
    additivity <- fina_result[fina_result$`abs(d/a)`<=0.2,]
    additivity$condition <- "Additivity"
  }else{
    additivity <- data.frame()
  }
  if (length(da[(da>0.2)&(da<=0.2)]>0)){
    partialDom <- fina_result[(fina_result$`abs(d/a)`>0.2)&(fina_result$`abs(d/a)`<=0.8),]
    partialDom$condition <- "Partial Dominance"
  }else{
    partialDom <- data.frame()
  }
  if (length(da[(da>0.8)&(da<=1.2)]>0)){
    Dominance <- fina_result[(fina_result$`abs(d/a)`>0.8)&(fina_result$`abs(d/a)`<=1.2),]
    Dominance$condition <- "Dominance"
  }else{
    Dominance <- data.frame()
  }
  if (length(da[da>1.2]>0)){
    Overdominance <- fina_result[fina_result$`abs(d/a)`>1.2,]
    Overdominance$condition <- "Overdominance"
  }else{
    Overdominance <- data.frame()
  }
  result <- rbind(additivity,partialDom,Dominance,Overdominance)
  row.names(result)<- 1:nrow(result)
  return(result)
}
