#structToolbox - https://github.com/computational-metabolomics/structToolbox/blob/0aafdd24cf4358e54610ce9b94001c7bcfdee72c/R/ggplot_theme_pub.R

#' Ggplot2 theme for publication figures
#'
#' @param base_size Set base size
#'
#' @return A ggplot2 theme object to add to a plot.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_Publication()
#'
#' @export theme_Publication
theme_Publication <- function(base_size=14){ #, base_family="helvetica") {
  (ggthemes::theme_foundation(base_size=base_size) #, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = '#ffffff'),
           plot.background = element_rect(colour = '#ffffff'),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = '#ffffff'),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.75, "cm"),
           legend.spacing = unit(0.2, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))

}


#' Ggplot2 pal for publication figures
#'
#' @param ... Arguments passed on to [ggplot2::discrete_scale()].
#'
#' @return A ggplot2 discrete fill scale.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(factor(cyl), fill = factor(cyl))) +
#'   geom_bar() + scale_fill_Publication()
#'
#' @export scale_fill_Publication
scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication",scales::manual_pal(values = c("#386cb0","#ef3b2c","#7fc97f","#fdb462","#984ea3","#a6cee3","#778899","#fb9a99","#ffff33")), ...)
}

#' Ggplot2 pal for publication figures
#'
#' @param ... Arguments passed on to [ggplot2::discrete_scale()].
#'
#' @return A ggplot2 discrete colour scale.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(factor(cyl), mpg, colour = factor(cyl))) +
#'   geom_point() + scale_colour_Publication()
#'
#' @export scale_colour_Publication
scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication",scales::manual_pal(values = c("#386cb0","#ef3b2c","#7fc97f","#fdb462","#984ea3","#a6cee3","#778899","#fb9a99","#ffff33")), ...)

}

#' Build an ordered class factor and matching colour palette
#'
#' @param class Vector of class labels (character or factor).
#' @param QC_label,Blank_label Labels treated as QC / Blank and pinned to the
#'   front of the ordering with fixed colours.
#' @param QC_color,Blank_color Colours for the QC / Blank classes.
#' @param manual_color Colour vector used for the remaining classes.
#'
#' @return A list with `class` (an ordered factor) and `manual_colors`.
#'
#' @examples
#' createClassAndColors(class = c("Sample", "Sample", "QC", "Blank"))
#'
#' @export createClassAndColors
createClassAndColors <- function (class, QC_label="QC", Blank_label="Blank", QC_color="#000000",
                                  Blank_color="#A65628",
                                  manual_color=c("#386cb0","#ef3b2c","#7fc97f","#fdb462","#984ea3","#a6cee3","#778899","#fb9a99","#ffff33")
)
{
  if (!is.ordered(class))
  {
    reorderNames <- sort(as.character(unique(class)))
  } else
  {
    reorderNames=levels(class)
  }

  # if too many levels then use rainbow scale with suitable number of colours
  if (length(reorderNames)>length(manual_color))
  {
    manual_color=grDevices::rainbow(length(reorderNames))
  }

  hit1 <- which(reorderNames==QC_label)
  if (length(hit1)==0) hit1 <- NULL

  hit2 <- which(reorderNames==Blank_label)
  if (length(hit2)==0) hit2 <- NULL

  if (!is.null(c(hit1,hit2)))
  {
    remo <- seq_along(reorderNames)[-c(hit1,hit2)]
  } else
  {
    remo <- seq_along(reorderNames)
  }

  reorderNames <- reorderNames[c(hit1, hit2,remo)]

  class <- factor (class, levels=reorderNames, ordered=TRUE)

  extraColors <- NULL
  if(!is.null(hit1)) extraColors[1] <- QC_color
  if(!is.null(hit2)) extraColors[2] <- Blank_color

  if (!is.null(extraColors)) extraColors <- extraColors[!is.na(extraColors)]

  out <- list(class=class,manual_colors=c(extraColors,manual_color))

  return(out)
}
