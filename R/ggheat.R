#' ggplot2 heatmap
#'
#' Draws a (potentially clustered) heatmap, with a lighter-weight resulting
#' figure file than is typical for heatmap.2.
#'
#' @param m A numeric matrix or data.frame
#' @param rescaling Rescale? Options are "row" or "column"
#' @param clustering Cluster entries? Options are "row", "column", or "both".
#' @param labCol Label columns?
#' @param labRow Label rows?
#' @param border Add border to plot?
#' @param heatscale Color range for colorscale.
#'
#' @return A ggplot object containing the formatted heatmap
#'
#' @details
#'
#'
#' @examples
#'
#' data(mtcars)
#' x=as.matrix(mtcars)
#' ggheat(x, clustering='column', rescaling='row', heatscale=c(low='red',
#' high='yellow'))
#'
#' @seealso \code{\link{runDiagnostics}}, \code{\link{plotCorrelations}}
#'
#' @export

ggheat=function(m, rescaling='none', clustering='none', labCol=T, labRow=T,
                border=FALSE, heatscale= c(low='darkblue',high='lightblue'))
{

  # This code adapted from a snippet found at:
  # https://www.r-bloggers.com/ggheat-a-ggplot2-style-heatmap-function/

  if(is.function(rescaling))
  {
    m=rescaling(m)
  }
  else
  {
    if(rescaling=='column')
      m=scale(m, center=T)
    if(rescaling=='row')
      m=t(scale(t(m),center=T))
  }

  if(is.function(clustering))
  {
    m=clustering(m)
  }else
  {
    if(clustering=='row')
      m=m[hclust(dist(m))$order, ]
    if(clustering=='column')
      m=m[,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }

  rows=dim(m)[1]
  cols=dim(m)[2]
  melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows),
               reshape::melt.array(m))
  g=ggplot2::ggplot(data=melt.m)


  if(border==TRUE)
    g2=g+ggplot2::geom_rect(ggplot2::aes(xmin=colInd-1,xmax=colInd,
                                         ymin=rowInd-1,ymax=rowInd, fill=value)
                            ,colour='white')
  if(border==FALSE)
    g2=g+ggplot2::geom_rect(ggplot2::aes(xmin=colInd-1,xmax=colInd,
                                         ymin=rowInd-1,ymax=rowInd, fill=value))

  if(labCol==T)
    g2=g2+ggplot2::scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
  if(labCol==F)
    g2=g2+ggplot2::scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))

  if(labRow==T)
    g2=g2+ggplot2::scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))
  if(labRow==F)
    g2=g2+ggplot2::scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))

  g2=g2+ggplot2::theme(panel.grid.minor=ggplot2::element_line(colour=NA),
                       panel.grid.major=ggplot2::element_line(colour=NA),
                       panel.background=ggplot2::element_rect(fill=NA,
                                                              colour=NA))

  return(g2+ggplot2::scale_fill_continuous(low=heatscale[1], high=heatscale[2],
                                           guide="colourbar"))

}

