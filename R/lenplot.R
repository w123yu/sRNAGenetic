#' Generate the sRNA length distribution plot
#'
#' @param file_dataframe A dataframe. The output result after running "srnapredata".
#' @param width A numeric. The width of the output bar plot, and default is 0.6.
#' @param size A numeric. The size of text in the outpot plot, and default is 12.
#'
#' @return The sRNA length distribution plot
#' @export
#' @examples
#' ##F1
#' F1_sRNA <- srnapredata(srnaseq_dataframe = F1_sRNA_seq, group = "F1")
#' ##P1
#' P1_sRNA <- srnapredata(srnaseq_dataframe = P1_sRNA_seq, group = "P1")
#' ##P2
#' P2_sRNA <- srnapredata(srnaseq_dataframe = P2_sRNA_seq, group = "P2")
#' ##integrate all sRNA data from P1, P2, and F1
#' sRNA_data <- rbind(F1_sRNA,P1_sRNA,P2_sRNA)
#' ##plot
#' lenplot(file_dataframe = sRNA_data)
lenplot <- function(file_dataframe,width = 0.6,size = 12){
  colnames(file_dataframe) <- c("Length","Frequency","Group")
  srna <- plyr::ddply(file_dataframe,"Group",transform,Percent = Frequency/sum(Frequency)*100)
  srna_plot <- ggplot2::ggplot(srna,ggplot2::aes(Length,Percent,fill=Group))+
    ggplot2::geom_bar(stat = "identity",position = "dodge",width = width)+
    ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
    ggplot2::theme(axis.text = ggplot2::element_text(size = size,family="serif"),
                   axis.title = ggplot2::element_text(size = size,face = "bold",family="serif"),
                   legend.text = ggplot2::element_text(size = size,family = "serif"),
                   legend.title = ggplot2::element_text(size = size,family = "serif",face = "bold"),
                   panel.border = ggplot2::element_rect(fill=NA,color="black", size=1, linetype="solid")
    )
  return(srna_plot)
}
utils::globalVariables(
  c("Length","Frequency","Group","Percent","Base","Position","colors","na.omit","write.csv")
)
