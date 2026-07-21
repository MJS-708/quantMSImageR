setGeneric("imageR", function(MSIobject, ...) standardGeneric("imageR"))

#' Draw an ion image for one feature
#'
#' Renders a single feature's spatial distribution as a `ggplot`. Because the
#' result is an ordinary `ggplot`, the layout, theme and scales can be modified
#' with the usual `+` syntax -- for example `+ facet_wrap(~ sample, ncol = 3)`
#' to arrange several acquisitions in a grid.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @param val_slot Character. Spectra slot to image (default `"intensity"`).
#' @param value Character. Legend label describing what the values represent
#'   (e.g. `"pg/mm2"` for a calibrated layer).
#' @param scale Character. Colour-scale treatment: `"suppress"` (cap at
#'   `percentile`), `"histogram"` (contrast-enhance per image) or `"sqrt"`.
#' @param percentile Numeric. Percentile at which the colour scale is capped
#'   when `scale = "suppress"` (default `99`).
#' @param threshold Numeric. Lowest percentile removed from the colour scale,
#'   to stop residual background dominating it (default `1`).
#' @param sample_lab Character. Column of `pData(MSIobject)` used to label and
#'   facet the images (default `"sample_ID"`).
#' @param pixels Character. Value of `pData(MSIobject)$sample_type` selecting
#'   which pixels to plot; `NA` (the default) plots all.
#' @param feat_ind Integer. Row index of the feature in `fData(MSIobject)` to
#'   image (default `1`).
#' @param perc_scale Logical. Express the colour scale as a percentage of the
#'   maximum (`TRUE`) rather than raw values (`FALSE`, the default).
#' @param blank_back Logical; when TRUE background/zero pixels are drawn transparent.
#' @param aspect_ratio numeric plot aspect ratio (default 1).
#' @param text_image Logical; when TRUE return the image as a numeric matrix
#'   (for text export) rather than a ggplot.
#' @return A `ggplot` object, or a numeric matrix when `text_image = TRUE`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' gg <- imageR(obj, feat_ind = 1, sample_lab = "run", scale = "suppress")
#'
#' @family visualisation
#' @aliases imageR
#' @export
setMethod("imageR", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", value = "response %", scale = "suppress", threshold = 1,
                   sample_lab = "sample_ID", pixels = NA, percentile=99.0,
                   feat_ind = 1, perc_scale = FALSE, blank_back = TRUE, aspect_ratio=1,
                   text_image = FALSE){

            MSIobject = as(MSIobject[feat_ind, ], "quant_MSImagingExperiment")

            if(!is.na(pixels)){
              MSIobject = MSIobject[, which(pData(MSIobject)$sample_type == pixels)]
            }

            # Generate image data matrix
            image_df = tibble::tibble(x = pData(MSIobject)@listData[["x"]],
                              y = pData(MSIobject)@listData[["y"]],
                              response = as.numeric(spectraData(MSIobject)[[val_slot]]),
                              sample = pData(MSIobject)[[sample_lab]],
                              feature = fData(MSIobject)$name)


            if(scale == "sqrt"){
              image_df = mutate(image_df, response = sqrt(response))
            }
            if(scale == "suppress"){
              vals = image_df$response
              vals[which(vals < 0)] = 0
              max_val = quantile(vals, percentile /100, na.rm=TRUE)
              if (!is.na(max_val) && max_val > min(vals, na.rm=TRUE)){
                vals[vals > max_val] = max_val
              }

              image_df$response = vals
            }
            if(scale == "histogram"){

              # Cardinal implementation https://rdrr.io/bioc/Cardinal/src/R/DIP.R

              vals = image_df$response

              unique_percentiles = unique(quantile(vals, seq(from=0, to=1, length.out=100), na.rm=TRUE))
              if (length(unique_percentiles) < 5){
                return("Range of intensity values too narrow for histogram normalization.") }

              vals_interval = cut(vals, unique_percentiles, include.lowest=TRUE)
              vals_new = as.numeric(vals_interval) / length(levels(vals_interval))

              scale = mean(vals, na.rm=TRUE) / mean(vals_new, na.rm=TRUE) # So mean value same as prior to scaling
              vals_new = scale * vals_new

              image_df$response = vals_new

            }

            if(!is.na(threshold)){
              # Remove x% of values - threshold
              vals = image_df$response
              min_val = quantile(vals, threshold /100, na.rm=TRUE)
              if (!is.na(min_val) && min_val < max(vals, na.rm=TRUE)){
                vals[vals < min_val] = 0
              }
              image_df$response = vals
            }

            # Set NA and negative values to 0 for plotting!
            image_df = image_df |>
              dplyr::mutate(response = ifelse(is.na(response), 0, response)) |>
              dplyr::mutate(response = ifelse(response < 0, 0, response))

            if(perc_scale == TRUE){
              vals = image_df$response
              if(max(vals) > 0){
                perc_vals = 100 * (vals / max(vals))

                image_df$response = perc_vals

              }
            }

            if(text_image == TRUE){

              p = image_df |>
                dplyr::select(x, y, response) |>
                dplyr::arrange(as.numeric(x)) |>
                tidyr::pivot_wider(names_from = x, values_from = response) |>
                dplyr::arrange(as.numeric(y)) |>
                tibble::column_to_rownames("y")

            } else if(blank_back == TRUE){
              image_df$response[is.na(image_df$response) | image_df$response <= 0] <- NA

              p = ggplot(data=image_df, aes(x = x, y = -y, fill = response)) +
                geom_tile() +
                theme_minimal() +
                theme(
                  aspect.ratio = aspect_ratio,
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.line = element_blank(),
                  panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
                ) +
                scale_fill_viridis(na.value = "white") +
                labs(fill = value) +
                facet_grid(sample ~ feature)


            } else{
              image_df[is.na(image_df)] <- 0

              p = ggplot(data=image_df,aes(x=x,y=-y,fill=response))+
                geom_tile() +
                theme_minimal() +
                theme(aspect.ratio=aspect_ratio,
                      axis.title = element_blank(),
                      axis.text = element_blank(),
                      axis.line = element_blank(),
                      panel.grid = element_blank(),
                      plot.title = element_text(hjust = 0.5, face="bold", size = 15)) +
                scale_fill_viridis(na.value = "white") +
                labs(fill=value) +
                facet_grid(sample~feature)
            }

            return(p)

        })
