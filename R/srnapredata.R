#' Generate the data of sRNA length distribution
#'
#' Generally, the length interval of sRNA is 21-24. The function of "srnapredata" can provide the input data for the next drawing of sRNA length distribution among different species.
#'
#' @param srnaseq_dataframe A dataframe. The first column must be the sRNA sequence.
#' @param group A character. You an select a representative group name for next drawing.
#'
#' @return A dataframe. The output results are consist of three columns, the first column is the length of sRNA, the second column id the frequency, and the third column is the group name.
#' @export
#' @examples
#' ##Only 400 sRNAs are selected as test data due to the large data of sRNA.
#' ##Recommended to use the "data.table" package for reading data quickly.
#' ##F1
#' F1_sRNA <- srnapredata(srnaseq_dataframe = F1_sRNA_seq, group = "F1")
#' ##P1
#' P1_sRNA <- srnapredata(srnaseq_dataframe = P1_sRNA_seq, group = "P1")
#' ##P2
#' P2_sRNA <- srnapredata(srnaseq_dataframe = P2_sRNA_seq, group = "P2")
srnapredata <- function(srnaseq_dataframe,group){
  srnaseq_dataframe$length <- nchar(as.vector(as.matrix(srnaseq_dataframe[,1])))
  res <- as.data.frame(table(subset(srnaseq_dataframe,select=c(length))))
  res$group <- group
  colnames(res) <- c("Length","Frequency","Group")
  return(res)
}
