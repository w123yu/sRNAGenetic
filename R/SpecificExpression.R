#' Species specific expression analysis: miVennPlot
#'
#' miVennPlot: generate the Venn diagram with the specific expression information of miRNAs.
#'
#' @param P1_RPM A dataframe. The rpm data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the rpm of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the P2 species.
#' @param F1_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the F1 species.
#' @param rpm_threshold A numeric. the average of rpm value among all the biological replicates. By default, the average of rpm more than or equal to 1 is retained.
#'
#' @return The Venn diagram with the specific expression information of miRNAs.
#' @export
#' @examples
#' ##Drawing
#' miVennPlot(P1_RPM = P1_miRNA_rpm,
#'            P2_RPM = P2_miRNA_rpm,
#'            F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
miVennPlot <- function(P1_RPM,P2_RPM,F1_RPM,rpm_threshold = 1){
  func_mirnafiliter <- function(data1,threshold){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- as.vector(as.matrix(data2[,1]))
    return(result)
  }
  P1_mirna <- func_mirnafiliter(data1 = P1_RPM, threshold = rpm_threshold)
  P2_mirna <- func_mirnafiliter(data1 = P2_RPM, threshold = rpm_threshold)
  F1_mirna <- func_mirnafiliter(data1 = F1_RPM, threshold = rpm_threshold)
  x = list(P1=P1_mirna,P2=P2_mirna,F1=F1_mirna)
  venn.plot <- VennDiagram::venn.diagram(
    x,euler.d = TRUE,filename = NULL,fontfamily = "serif",col="white",
    fill=c(colors()[616], colors()[38], colors()[468]),
    alpha=c(0.4, 0.4, 0.4),lwd=c(0.1, 0.1, 0.1),
    cat.dist = c(0.03, 0.03, 0.03), #调整分类条目位置
    cex = 2,cat.cex = 2,reverse = TRUE);
  grid::grid.draw(venn.plot)
}

#' Species specific expression analysis: miVennData
#'
#' miVennData: Extract the species-specific miRNAs and the shared miRNAs among parents and offspring.
#'
#' @param P1_RPM A dataframe. The rpm data of miRNA from the P1 species. The first column must be the miRNA sequence. Others are listed as the rpm of miRNA, and each column denotes one biological replicate of the sample.
#' @param P2_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the P2 species.
#' @param F1_RPM A dataframe. Similar with P1_RPM, the rpm data of miRNA from the F1 species.
#' @param rpm_threshold A numeric. the average of rpm value among all the biological replicates. By default, the average of rpm more than or equal to 1 is retained.
#' @param output_file Specify the output file. "venn_list" is the default option, which outputs all the information of the Venn diagram. "all_common" is one of options, which outputs the miRNAs shared by parents and offspring. "P1_specific" is one of options, which outputs P1 specific expression miRNA. "P2_specific" is one of options, which outputs P2 specific expression miRNA. "F1_specific" is one of options, which outputs F1 specific expression miRNA.
#'
#' @return A dataframe. The output results is based on your selection (output_file).
#' @export
#'
#' @examples
#' ##Extract the species-specific miRNAs and the shared miRNAs among parents and offspring.
#' ##output_file = "venn_list"
#' venn_list <- miVennData(P1_RPM = P1_miRNA_rpm,
#'                         P2_RPM = P2_miRNA_rpm,
#'                         F1_RPM = F1_miRNA_rpm,
#'                         rpm_threshold = 1,output_file = "venn_list")
#' ##output_file = "P1_specific"
#' P1_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           rpm_threshold = 1,output_file = "P1_specific")
#' ##output_file = "P2_specific"
#' P2_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           rpm_threshold = 1,output_file = "P2_specific")
#' ##output_file = "F1_specific"
#' F1_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
#'                           P2_RPM = P2_miRNA_rpm,
#'                           F1_RPM = F1_miRNA_rpm,
#'                           rpm_threshold = 1,output_file = "F1_specific")
#' ##output_file = "all_common"
#' all_common <- miVennData(P1_RPM = P1_miRNA_rpm,
#'                          P2_RPM = P2_miRNA_rpm,
#'                          F1_RPM = F1_miRNA_rpm,
#'                          rpm_threshold = 1,output_file = "all_common")
miVennData <- function(P1_RPM,P2_RPM,F1_RPM,rpm_threshold = 1,output_file="venn_list"){
  func_mirnafiliter <- function(data1,threshold){
    Average_rpm <- as.data.frame(apply(data1[,c(-1)],1,mean))
    data1_value <- Average_rpm >= threshold
    data2 <- data1[data1_value,]
    result <- as.vector(as.matrix(data2[,1]))
    return(result)
  }
  P1_mirna <- func_mirnafiliter(data1 = P1_RPM,threshold = rpm_threshold)
  P2_mirna <- func_mirnafiliter(data1 = P2_RPM,threshold = rpm_threshold)
  F1_mirna <- func_mirnafiliter(data1 = F1_RPM,threshold = rpm_threshold)
  x = list(P1=P1_mirna,P2=P2_mirna,F1=F1_mirna)
  inter <- VennDiagram::get.venn.partitions(x)
  for (i in 1:nrow(inter)){
    inter[i,'..values..'] <- paste(inter[[i,'..values..']], collapse = ',')
  }
  for (i in 1:nrow(inter)){
    if (inter[i,1]=="TRUE" & inter[i,2]=="TRUE" & inter[i,3]=="TRUE"){
      P1_P2_F1 <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P1_P2_F1) <- "P1_P2_F1"
    }else if (inter[i,1]=="TRUE" & inter[i,2]=="FALSE" & inter[i,3]=="FALSE"){
      P1_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P1_only) <- "P1_only"
    }else if (inter[i,2]=="TRUE" & inter[i,1]=="FALSE" & inter[i,3]=="FALSE"){
      P2_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(P2_only) <- "P2_only"
    }else if (inter[i,3]=="TRUE" & inter[i,1]=="FALSE" & inter[i,2]=="FALSE"){
      F1_only <- as.data.frame(strsplit(as.character(inter[i,'..values..']),","))
      names(F1_only) <- "F1_only"
    }else{
      next
    }
  }
  if (output_file == "all_common"){
    return(P1_P2_F1)
  }else if (output_file == "P1_specific"){
    return(P1_only)
  }else if (output_file == "P2_specific"){
    return(P2_only)
  }else if (output_file == "F1_specific"){
    return(F1_only)
  }else{
    inter$..values.. <- unlist(inter$..values..)
    colnames(inter) <- c("P1","P2","F1","set","values","count")
    return(inter)
  }
}
