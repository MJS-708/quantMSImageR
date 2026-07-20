setGeneric("createMSIDatamatrix", function(MSIobject, ...) standardGeneric("createMSIDatamatrix"))

#' Function to create data matrix from MSI object
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSI object from Cardinal, pData() to include sample_ID.
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param inputNA Whether to convert 0's in matrix to NA (default = TRUE)
#' @param roi_header Header in pData pertaining to ROIs to average. Set to NA to skip generating average df
#' @return MSIobject with slots updated for i) matrix of average ng/pixel of m/z (rows = m/z and cols = cal level) in tissue ROIs ii) sample/ROI metadata
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- createMSIDatamatrix(obj, val_slot = "intensity", roi_header = NA)
#'
#' @aliases createMSIDatamatrix
#' @export
setMethod("createMSIDatamatrix", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", inputNA = TRUE, roi_header = NA){

            # Subset pixels in ROIs only (based on roi_header). ROI is the
            # grouping key used below: the pixel index when no roi_header is
            # given, otherwise the roi_header column itself.
            if(is.na(roi_header)){
              pData(MSIobject)$ROI = seq_len(nrow(pData(MSIobject)))
            } else{
              MSIobject = MSIobject[, which(!is.na(pData(MSIobject)[[roi_header]]))]
              pData(MSIobject)$ROI = pData(MSIobject)[[roi_header]]
            }

            # Update pixel data
            pixel_df = data.frame(pData(MSIobject)) |>
              subset(!is.na(ROI)) |>
              tibble::rownames_to_column("pixel_ind") |>
              dplyr::mutate(pixel_ind = sprintf("pixel_%s", pixel_ind))


            # All pixel df
            all_pixel_df = data.frame(
              vapply(seq_len(nrow(fData(MSIobject))),
                     function(x) spectraData(MSIobject)[[val_slot]][x, ],
                     numeric(ncol(MSIobject)))
            ) |>
              dplyr::mutate(pixel_ind = pixel_df$pixel_ind)
            colnames(all_pixel_df) = c(fData(MSIobject)$name, "pixel_ind")

            .keep_cols = c(fData(MSIobject)$name, "pixel_ind",
                           if (!is.na(roi_header)) roi_header)
            all_pixel_df = dplyr::left_join(all_pixel_df, pixel_df, by = "pixel_ind") |> dplyr::select(dplyr::any_of(.keep_cols))


            if(inputNA){
              all_pixel_df <- replace(all_pixel_df, all_pixel_df==0, NA)
            }


            # Create average ROI df
            if(!is.na(roi_header)){

              #### THIS NEEDS FIXING TO SUMMARISE EACH FEATURE INDEPENDENTLY!!!!
              ave_df = all_pixel_df |>
                dplyr::group_by(dplyr::across(dplyr::all_of(roi_header))) |>
                dplyr::summarise(dplyr::across(dplyr::any_of(c(fData(MSIobject)$name)), \(x) mean(x, na.rm = TRUE))) |>
                tibble::column_to_rownames(roi_header)

              if(inputNA){
                ave_df <- replace(ave_df, ave_df==0, NA)
              }

              MSIobject@tissueInfo@roi_average_matrix = ave_df
            }

            all_pixel_df = all_pixel_df |>
              tibble::column_to_rownames("pixel_ind") |>
              dplyr::select(dplyr::any_of(c(fData(MSIobject)$name)))

            MSIobject@tissueInfo@all_pixel_matrix = all_pixel_df
            MSIobject@tissueInfo@sample_metadata = pixel_df

            return(MSIobject)
          })
