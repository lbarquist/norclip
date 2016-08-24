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

clip_scat <- function(exp, ctrl, sf_exp=1, sf_ctrl=1, xyline=F, exp_lab="crosslinked", ctrl_lab="non-crosslinked", main="2D density plot", elliptical=2, nbin=100, colramp=rev(rainbow(10, end = 4/6))){
  if(elliptical > 0){
    filt <- filter_elliptical(exp, ctrl, elliptical)
    masked_exp <- as.vector(exp[filt], mode="integer")
    masked_ctrl <- as.vector(ctrl[filt], mode="integer")
  }else{
    masked_exp <- as.vector(exp, mode="integer")
    masked_ctrl <- as.vector(ctrl, mode="integer")
  }
  masked_exp <- masked_exp / sf_exp
  masked_ctrl <- masked_ctrl / sf_ctrl

  this_dat <- data.frame(exp = masked_exp, ctrl = masked_ctrl)

  p <- ggplot(this_dat, aes(x=exp, y=ctrl))
  p <- p + geom_hex(aes(x=exp,y=ctrl), bins=100) + xlab("+XL") + ylab ("-XL")
  p <- p +scale_fill_gradientn("", colours = colramp) + theme_bw()
  p <- p + ggtitle(main)
  if(xyline){
    p <- p + geom_abline(intercept=0, slope=1, colour="red")
  }
  print(p)
  #smoothScatter(masked_exp, masked_ctrl, nbin=nbin, nrpoints=nrpoints,colramp=colramp, pch=19, cex=.4, xlab=exp_lab, ylab=ctrl_lab)
}
