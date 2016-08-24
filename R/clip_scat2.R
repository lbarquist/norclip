#' A function for producing smoothed scatterplots for CLIP-seq data
#'
#'
#' @param exp A vector containing read-counts per nucleotide in a crosslinked sample
#' @param ctrl A vector containing read-counts per nucleotide in a non-crosslinked control
#' @param exp_lab A label describing the experimental condition
#' @param ctrl_lab A label describing the control condition
#' @param elliptical The number of standard deviations to use in constructing an elliptical filter
#'
#' @examples
#'
#' @export
#'
#'

clip_scat2 <- function(exp, ctrl, exp_lab="crosslinked", ctrl_lab="non-crosslinked", elliptical=2, nbin=1000, nrpoints=0, colramp=colorRampPalette(rev(brewer.pal(11, "RdYlBu")))){
  filt <- filter_elliptical(exp, ctrl, elliptical)
  masked_exp <- as.vector(exp[filt])
  masked_ctrl <- as.vector(ctrl[filt])
  df <- data.frame(exp=masked_exp, ctrl=masked_ctrl)
  ggplot(df,aes(x=x,y=y))+
    stat_density2d(aes(alpha=..level..), geom="polygon") +
    scale_alpha_continuous(limits=c(0,0.2),breaks=seq(0,0.2,by=0.025))+
    geom_point(colour="red",alpha=0.02)+
    theme_bw()
}
