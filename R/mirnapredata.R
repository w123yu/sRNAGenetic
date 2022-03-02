#' Generate the data of miRNA base frequence in each position
#'
#' Generally, the "T" base account for the highest percentage of miRNA in the first position.The function of "mirnapredata" can provide the input data for the next drawing of miRNA base distribution in each position.
#'
#' @param mirnaseq_dataframe A dataframe. The first column must be the sRNA sequence.
#'
#' @return A dataframe. About the output results, the first column is the base, the second column is the base frequency, the third column is the position.
#' @export
#'
#' @examples
#' ##P1
#' P1_miRNA_data <- mirnapredata(mirnaseq_dataframe = P1_miRNA_count)
#' ##P2
#' P2_miRNA_data <- mirnapredata(mirnaseq_dataframe = P2_miRNA_count)
#' ##F1
#' F1_miRNA_data <- mirnapredata(mirnaseq_dataframe = F1_miRNA_count)
mirnapredata <- function(mirnaseq_dataframe){
  max_len <- max(nchar(as.vector(as.matrix(mirnaseq_dataframe[,1]))))
  mirna_list <- strsplit(as.vector(as.matrix(mirnaseq_dataframe[,1])),"")
  new_list <- lapply(mirna_list, function(x) {c(x, rep(NA, max_len - length(x)))})
  mirna_matrix <- do.call(rbind,new_list)
  res <- data.frame()
  for(i in (1:max_len)){
    data <- as.data.frame(table(mirna_matrix[,i]))
    data$Position <- i
    res <- rbind(res,data)
  }
  names(res) <- c("Base","Frequency","Position")
  return(res)
}
