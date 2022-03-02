#' Generate the base frequency plot of miRNA
#'
#' @param file_dataframe A dataframe. The output result after running mirnapredata.
#' @param width A numeric. The width of the output bar plot, and default is 0.6.
#' @param size A numeric. The size of axis text, and default is 0.6.
#'
#' @return The miRNA base frequency plot
#' @export
#'
#' @examples
#' ##P1
#' P1_miRNA_data <- mirnapredata(mirnaseq_dataframe = P1_miRNA_count)
#' ##P2
#' P2_miRNA_data <- mirnapredata(mirnaseq_dataframe = P2_miRNA_count)
#' ##F1
#' F1_miRNA_data <- mirnapredata(mirnaseq_dataframe = F1_miRNA_count)
#' ##Drawing
#' basepreplot(file_dataframe = P1_miRNA_data)
#' basepreplot(file_dataframe = P2_miRNA_data)
#' basepreplot(file_dataframe = F1_miRNA_data)
basepreplot <- function(file_dataframe,width = 0.6,size = 12){
  mirna <- plyr::ddply(file_dataframe,"Position",transform,Percent = Frequency/sum(Frequency)*100)
  ggplot2::ggplot(mirna,ggplot2::aes(Position,Percent,fill=Base))+
    ggplot2::scale_x_continuous(breaks = seq(min(mirna$Position),max(mirna$Position),1))+
    ggplot2::geom_bar(stat = "identity",position = "stack",width = width)+
    ggsci::scale_fill_npg()+ggplot2::theme_classic()+ggplot2::xlab("Position")+ggplot2::ylab("Percent")+
    ggplot2::theme(axis.text = ggplot2::element_text(size = size,family="serif"),
                   axis.title = ggplot2::element_text(size = size,face = "bold",family="serif"),
                   legend.text = ggplot2::element_text(size = size,family = "serif"),
                   legend.title = ggplot2::element_text(size = size,family = "serif",face = "bold"),
                   panel.border = ggplot2::element_rect(fill=NA,color="black", size=1, linetype="solid")
    )
}
