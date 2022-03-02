#' Classification for 12 expression patterns
#'
#' The input data is generated from the analysis result of DESeq2.
#'
#' @param pv11 A numeric. The P value of F1_vs_P1 (Treatment:F1; Control:P1).
#' @param pv12 A numeric. The P value of F1_vs_P2 (Treatment:F1; Control:P2)
#' @param pv21 A numeric. The P value of P2_vs_P1 (Treatment:P2; Control:P1)
#' @param fc11 A numeric. The Log2(FoldChange) value of F1_vs_P1 (Treatment:F1; Control:P1)
#' @param fc12 A numeric. The Log2(FoldChange) value of F1_vs_P2 (Treatment:F1; Control:P2)
#' @param fc21 A numeric. The Log2(FoldChange) value of P2_vs_P1 (Treatment:P2; Control:P1)
#' @param Pvalue A numeric. Filtration criteria of P value for Classification.
#'
#' @return A dataframe.
#' @export
#'
genetic <- function(pv11,pv12,pv21,fc11,fc12,fc21,Pvalue){
  if(pv11 < Pvalue && pv12 < Pvalue){
    if(fc11 < -1 && fc12 > 1 && pv21 < Pvalue){
      return("XII")
    }else if(fc11 >= 1 && fc12 <= -1 && pv21 < Pvalue){
      return("I")
    }else if(fc11 >= 1 && fc12 >= 1 && pv21 >= Pvalue){
      return("VIII")
    }else if(fc11 <= -1 && fc12 <= -1 && pv21 >= Pvalue){
      return("VII")
    }else if(fc11 >= 1 && fc12 >= 1 && pv21 < Pvalue){
      if(fc21 <= -1){
        return("VI")
      }else if(fc21 >= 1){
        return("V")
      }else{
        return("No change")
      }
    }else if(fc11 <= -1 && fc12 <= -1 && pv21 < Pvalue){
      if(fc21 <= -1){
        return("X")
      }else if (fc21 >= 1){
        return("III")
      }else{
        return("No change")
      }
    }
  }else if(pv11 >= Pvalue && pv12 < Pvalue && pv21 < Pvalue){
    if(fc21 <= -1 && fc12 >= 1){
      return("IV")
    }else if(fc21 >= 1 && fc12 <= -1){
      return("IX")
    }else{
      return("No change")
    }
  }else if(pv11 < Pvalue && pv12 >= Pvalue && pv21 < Pvalue){
    if(fc21 >= 1 && fc11 >= 1){
      return("II")
    }else if(fc21 <= -1 && fc11 <= -1){
      return("XI")
    }else{
      return("No change")
    }
  }else{
    return("No change")
  }
}
