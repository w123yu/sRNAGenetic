#' Filitering low expressed miRNAs based on count: Countfiliter
#'
#' @param P1_count A dataframe. The count data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the count of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_count A dataframe. Similar with P1_count, the count data of miRNA from the P2 species.
#' @param F1_count A dataframe. Similar with P1_count, the count data of miRNA from the F1 species.
#' @param count_threshold A numeric. In all samples, there is at least one sample whose count value is more than or equal to count_threshold to be retained. By default, the count value more than or equal to 5 is retained.
#'
#' @return A dataframe. The result includes all miRNAs that fulfill the count value requirement (count >= count_threshold) in at least one sample.
#' @export
#'
#' @examples
#' ##Get the filitered mirna count table (default: Count >= 5 in at least one sample)
#' Count5result <- Countfiliter(P1_count = P1_miRNA_count,
#'                              P2_count = P2_miRNA_count,
#'                              F1_count = F1_miRNA_count,count_threshold = 5)
Countfiliter <- function(P1_count,P2_count,F1_count,count_threshold = 5){
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
  rownames(filter_result) <- 1:nrow(filter_result)
  return(filter_result)
}

#' Filitering low expressed miRNAs based on RPM: Rpmfiliter
#'
#' @param P1_RPM A dataframe. The rpm data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the rpm of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the P2 species.
#' @param F1_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the F1 species.
#' @param rpm_threshold A numeric. the average of rpm value among all the biological replicates. By default, the average of rpm more than or equal to 1 is retained.
#'
#' @return A dataframe. The result includes all miRNAs that fulfill the average rpm value requirement (Average rpm >= rpm_threshold) among all species.
#' @export
#'
#' @examples
#' ##Get the filitered mirna rpm table (default: the average rpm >= 1 in three species)
#' Rpm1result <- Rpmfiliter(P1_RPM = P1_miRNA_rpm,
#'                          P2_RPM = P2_miRNA_rpm,
#'                          F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
Rpmfiliter <- function(P1_RPM,P2_RPM,F1_RPM,rpm_threshold = 1){
  if(ncol(P1_RPM) == ncol(P2_RPM) & ncol(P1_RPM) == ncol(F1_RPM)){
    colnum = ncol(P2_RPM)
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
  names(P1_RPM) <- P1_colname
  names(P2_RPM) <- P2_colname
  names(F1_RPM) <- F1_colname
  ############################
  func_mirnafiliter <- function(data){
    Average_rpm <- as.data.frame(apply(data[,c(-1)],1,mean))
    data1_value <- Average_rpm >= rpm_threshold
    data2 <- data[data1_value,]
    return(data2)
  }
  P1_mirna <- func_mirnafiliter(data = P1_RPM)
  P2_mirna <- func_mirnafiliter(data = P2_RPM)
  F1_mirna <- func_mirnafiliter(data = F1_RPM)
  ####################################
  Parent <- merge.data.frame(P1_mirna,P2_mirna,all=TRUE)
  filter_result <- merge.data.frame(Parent,F1_mirna,all=TRUE)
  fina_result <- na.omit(filter_result)
  ####################################
  rownames(fina_result) <- 1:nrow(fina_result)
  return(fina_result)
}
