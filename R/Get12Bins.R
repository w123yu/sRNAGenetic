#' Genetic effects analysis: Twelve bins of expression analysis (method2)
#'
#' @param P1_count A dataframe. The count data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the count of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_count A dataframe. Similar with P1_count, the count data of miRNA from the P2 species.
#' @param F1_count A dataframe. Similar with P1_count, the count data of miRNA from the F1 species.
#' @param count_threshold A numeric. In all samples, there is at least one sample whose count value is more than or equal to count_threshold to be retained. By default, the count value more than or equal to 5 is retained.
#' @param Pvalue A numeric. The threshold of significance test among different groups. Default is 0.05.
#'
#' @return A dataframe. The output results contain the P value, log2FoldChange and grouping information for each miRNA expressed in all species (count >= count_threshold). F1_vs_P1(P value: pv11,log2FoldChange: fc11), F1_vs_P2(P value: pv12,log2FoldChange: fc12), P2_vs_P1(P value: pv21,log2FoldChange: fc21)
#' @export
#'
#' @examples
#' ##Get the table of 12 expression patterns
#' Binresult <- Get12Bins(P1_count = P1_miRNA_count,
#'                        P2_count = P2_miRNA_count,
#'                        F1_count = F1_miRNA_count,
#'                        count_threshold = 5,Pvalue = 0.05)
Get12Bins <- function(P1_count,P2_count,F1_count,count_threshold = 5,Pvalue = 0.05){
  if(ncol(P1_count) == ncol(P2_count) & ncol(P1_count) == ncol(F1_count)){
    colnum = ncol(P2_count)
    sample_number = colnum-1
  }else{
    cat("Error!!! Inconsistent biological replicates between different samples.")
  }
  P1_colname <- "sequence"
  P2_colname <- "sequence"
  F1_colname <- "sequence"
  for (i in (1:sample_number)){
    P1_id <- paste("P1_",i,sep = "")
    P2_id <- paste("P2_",i,sep = "")
    F1_id <- paste("F1_",i,sep = "")
    P1_colname <- c(P1_colname,P1_id)
    P2_colname <- c(P2_colname,P2_id)
    F1_colname <- c(F1_colname,F1_id)
  }
  names(P1_count) <- P1_colname
  names(P2_count) <- P2_colname
  names(F1_count) <- F1_colname
  ############
  sum <- 3*sample_number
  func_countfiliter <- function(x_array){
    res <- "FALSE"
    for(i in 1:sum){
      if (x_array[i] >= count_threshold) {
        res <- "TRUE"
      } else {
        next
      }
    }
    if (res == "TRUE"){
      return("TRUE")
    }else{
      return("FALSE")
    }
  }
  input_data <- merge.data.frame(merge.data.frame(P1_count,P2_count,all=TRUE),F1_count,all=TRUE)
  input_data[is.na(input_data)] <- 0
  data1 <- input_data[,c(-1)]
  for(i in 1:nrow(data1)){
    if(i == 1){
      value_array <- func_countfiliter(as.vector(as.matrix(data1[1,])))
    }else{
      value_array <- c(value_array,func_countfiliter(as.vector(as.matrix(data1[i,]))))
    }
  }
  filter_result <- input_data[as.logical(value_array),]
  ####################################
  rownames(filter_result) <- filter_result[,1]
  filter_result <- filter_result[,c(-1)]
  ####################################
  P1.tab <- filter_result[,c(1:sample_number)]
  P2.tab <- filter_result[,c((sample_number+1):(2*sample_number))]
  F1.tab <- filter_result[,c((2*sample_number+1):(sample_number*3))]
  #############################
  P1_condition <- factor(c(rep("P1",sample_number),rep("F1",sample_number)),levels=c("P1","F1"))
  P2_condition <- factor(c(rep("P2",sample_number),rep("F1",sample_number)),levels=c("P2","F1"))
  P12_condition <- factor(c(rep("P1",sample_number),rep("P2",sample_number)),levels=c("P1","P2"))
  #############################
  P1F1.tab <- cbind(P1.tab,F1.tab)
  P2F1.tab <- cbind(P2.tab,F1.tab)
  P12.tab <- cbind(P1.tab,P2.tab)
  ##################################
  P1_colData <- data.frame(colnames(P1F1.tab),P1_condition)
  P2_colData <- data.frame(colnames(P2F1.tab),P2_condition)
  P12_colData <- data.frame(colnames(P12.tab),P12_condition)
  #############################
  P1_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P1F1.tab,colData = P1_colData,design = ~P1_condition)
  P2_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P2F1.tab,colData = P2_colData,design = ~P2_condition)
  P12_dds <- DESeq2::DESeqDataSetFromMatrix(countData = P12.tab,colData = P12_colData,design = ~P12_condition)
  ##############################
  cat("F1_vs_P1","\n")
  P1_dds <- DESeq2::DESeq(P1_dds)
  cat("F1_vs_P2","\n")
  P2_dds <- DESeq2::DESeq(P2_dds)
  cat("P2_vs_P1","\n")
  P12_dds <- DESeq2::DESeq(P12_dds)
  ########
  F1_vs_P1 <- as.data.frame(DESeq2::results(P1_dds))
  F1_vs_P2 <- as.data.frame(DESeq2::results(P2_dds))
  P2_vs_P1 <- as.data.frame(DESeq2::results(P12_dds))
  new_F1_vs_P1 <- cbind(as.data.frame(row.names(F1_vs_P1)),F1_vs_P1$pvalue,F1_vs_P1$log2FoldChange)
  new_F1_vs_P2 <- cbind(as.data.frame(row.names(F1_vs_P2)),F1_vs_P2$pvalue,F1_vs_P2$log2FoldChange)
  new_P2_vs_P1 <- cbind(as.data.frame(row.names(P2_vs_P1)),P2_vs_P1$pvalue,P2_vs_P1$log2FoldChange)
  names(new_F1_vs_P1) <- c("sequence","pv11","fc11")
  names(new_F1_vs_P2) <- c("sequence","pv12","fc12")
  names(new_P2_vs_P1) <- c("sequence","pv21","fc21")
  #####################
  SeGedata <- merge.data.frame(merge.data.frame(new_F1_vs_P1,new_F1_vs_P2,all=TRUE),new_P2_vs_P1,all=TRUE)
  SeGedata <- na.omit(SeGedata)
  va <- genetic(pv11=SeGedata[1,2],pv12=SeGedata[1,4],pv21=SeGedata[1,6],fc11=SeGedata[1,3],fc12=SeGedata[1,5],fc21=SeGedata[1,7],Pvalue = Pvalue)
  for(i in (2:nrow(SeGedata))){
    va2 <- genetic(pv11=SeGedata[i,2],pv12=SeGedata[i,4],pv21=SeGedata[i,6],fc11=SeGedata[i,3],fc12=SeGedata[i,5],fc21=SeGedata[i,7],Pvalue = Pvalue)
    va <- c(va,va2)
  }
  SeGedata$Categories <- va
  group_array <- levels(as.data.frame(table(SeGedata$Categories))$Var1)
  if("I" %in% group_array | "XII" %in% group_array){
    Additivity <- SeGedata[(SeGedata$Categories == "I") | (SeGedata$Categories == "XII"),]
    Additivity$Group <- "Additivity"
  }else{
    Additivity <- data.frame()
  }
  if("IV" %in% group_array | "IX" %in% group_array){
    ELD_P1 <- SeGedata[(SeGedata$Categories == "IV") | (SeGedata$Categories == "IX"),]
    ELD_P1$Group <- "P1-expression level dominance"
  }else{
    ELD_P1 <- data.frame()
  }
  if("II" %in% group_array | "XI" %in% group_array){
    ELD_P2 <- SeGedata[(SeGedata$Categories == "II") | (SeGedata$Categories == "XI"),]
    ELD_P2$Group <- "P2-expression level dominance"
  }else{
    ELD_P2 <- data.frame()
  }
  if("III" %in% group_array | "VII" %in% group_array | "X" %in% group_array ){
    TD <- SeGedata[(SeGedata$Categories == "III") | (SeGedata$Categories == "VII") | (SeGedata$Categories == "X"),]
    TD$Group <- "Transgressive downregulation"
  }else{
    TD <- data.frame()
  }
  if("V" %in% group_array | "VI" %in% group_array | "VIII" %in% group_array ){
    TU <- SeGedata[(SeGedata$Categories == "V") | (SeGedata$Categories == "VI") | (SeGedata$Categories == "VIII"),]
    TU$Group <- "Transgressive upregulation"
  }else{
    TU <- data.frame()
  }
  if("No change" %in% group_array){
    No <- SeGedata[SeGedata$Categories=="No change",]
    No$Group <- "No change"
  }else{
    No <- data.frame()
  }
  result <- rbind(Additivity,ELD_P2,ELD_P1,TD,TU,No)
  row.names(result) <- 1:nrow(result)
  return(result)
}
