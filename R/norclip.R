#' Perform norclip normalization of CLIP-seq libraries
#'
#'
#' @param file Path to a description file
#'
#' @examples
#'
#' @export
#'
#'

norclip <- function(file, crossnormalize=T, sdn=8, diagnostics=T, bg_cut=5, naive=T, export_wigs=F, wig_path="."){
  data_table <-read.table(file, stringsAsFactors = F)
  colnames(data_table) <- c("identifier","type","direction","file")
  wigs <- loadData(data_table)

  if(diagnostics){
    runDiagnostics(wigs, data_table, sdn=sdn)
  }

  factors <- clipScaleFactors(wigs, data_table, sdn, crossnormalize, plot=diagnostics, bg_cut=bg_cut)

  if(naive){
    naive <- naiveScaleFactors(wigs, data_table, sdn=sdn, bg_cut=bg_cut, plot=diagnostics)
  }

  return(factors)
}
