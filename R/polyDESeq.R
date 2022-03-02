#' @title Differential expression analysis
#'
#' @param P1_count A dataframe. The count data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the count of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_count A dataframe. Similar with P1_count, the count data of miRNA from the P2 species.
#' @param F1_count A dataframe. Similar with P1_count, the count data of miRNA from the F1 species.
#' @param count_threshold A numeric. In all samples, there is at least one sample whose count value is more than or equal to count_threshold to be retained. By default, the count value more than or equal to 5 is retained.
#' @param Pvalue A numeric. The threshold of significance test among different groups. Default is 0.05.
#'
#' @return A dataframe. Differential expression analysis results of miRNA expressed in each two species (count >= count_threshold).
#' @export
#'
#' @examples
#' ##Drawing
#' polyDESeq(P1_count = P1_miRNA_count,
#'           P2_count = P2_miRNA_count,
#'           F1_count = F1_miRNA_count,
#'           count_threshold = 5,Pvalue = 0.05)
polyDESeq <- function(P1_count,P2_count,F1_count,count_threshold = 5,Pvalue = 0.05){
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
  cat("DESeq results have been export to current directory:","\n","F1_vs_P1.csv,F1_vs_P2.csv,P2_vs_P1.csv","\n")
  ###########################
  F1_vs_P1 <- na.omit(F1_vs_P1)
  F1_vs_P2 <- na.omit(F1_vs_P2)
  P2_vs_P1 <- na.omit(P2_vs_P1)
  ####F1_vs_P1
  F1_vs_P1_up <-  F1_vs_P1[(F1_vs_P1$log2FoldChange > 1 & F1_vs_P1$pvalue < Pvalue), ]
  P1_vs_F1_up <-  F1_vs_P1[(F1_vs_P1$log2FoldChange < -1 & F1_vs_P1$pvalue < Pvalue), ]
  F1_vs_P1_num <- nrow(F1_vs_P1_up)
  P1_vs_F1_num <- nrow(P1_vs_F1_up)
  ####F1_vs_P2
  F1_vs_P2_up <-  F1_vs_P2[(F1_vs_P2$log2FoldChange > 1 & F1_vs_P2$pvalue < Pvalue), ]
  P2_vs_F1_up <-  F1_vs_P2[(F1_vs_P2$log2FoldChange < -1 & F1_vs_P2$pvalue < Pvalue), ]
  F1_vs_P2_num <- nrow(F1_vs_P2_up)
  P2_vs_F1_num <- nrow(P2_vs_F1_up)
  ####P1_vs_P2
  P2_vs_P1_up <-  P2_vs_P1[(P2_vs_P1$log2FoldChange > 1 & P2_vs_P1$pvalue < Pvalue), ]
  P1_vs_P2_up <-  P2_vs_P1[(P2_vs_P1$log2FoldChange < -1 & P2_vs_P1$pvalue < Pvalue), ]
  P2_vs_P1_num <- nrow(P2_vs_P1_up)
  P1_vs_P2_num <- nrow(P1_vs_P2_up)
  ###
  grid::grid.newpage()
  grid::grid.lines(x = grid::unit(c(1.5/5, 2.5/5), "npc"), y = grid::unit(c(1.5/5, 3.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL,gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.lines(x = grid::unit(c(1.5/5, 3.5/5), "npc"),y = grid::unit(c(1.5/5, 1.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL,gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.lines(x = grid::unit(c(3.5/5, 2.5/5), "npc"),y = grid::unit(c(1.5/5, 3.5/5), "npc"),
                   default.units = "npc",arrow = NULL, name = NULL, gp=grid::gpar(), draw = TRUE, vp = NULL)
  grid::grid.circle(x=1.5/5, y=1.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#0099FF"), draw=TRUE, vp=NULL)
  grid::grid.circle(x=3.5/5, y=1.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#FF9900"), draw=TRUE, vp=NULL)
  grid::grid.circle(x=2.5/5, y=3.5/5, r=0.1, default.units="npc", name=NULL, gp=grid::gpar(fill="#00CC00"), draw=TRUE, vp=NULL)
  grid::grid.text("P1",x=1.5/5,y=1.5/5)
  grid::grid.text("P2",x=3.5/5,y=1.5/5)
  grid::grid.text("F1",x=2.5/5,y=3.5/5)
  grid::grid.text(P1_vs_F1_num,x=1.65/5,y=2.15/5)
  grid::grid.text(P1_vs_P2_num,x=2/5,y=1.3/5)
  grid::grid.text(P2_vs_P1_num,x=3/5,y=1.3/5)
  grid::grid.text(P2_vs_F1_num,x=3.35/5,y=2.15/5)
  grid::grid.text(F1_vs_P1_num,x=2/5,y=3/5)
  grid::grid.text(F1_vs_P2_num,x=3/5,y=3/5)
}
