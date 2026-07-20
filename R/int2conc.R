setGeneric("int2conc", function(MSIobject, ...) standardGeneric("int2conc"))

#' Function to update intensity with concentration values
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject quant_MSImagingExperiment object
#' @param pixels Label in pixel metadata the pixels to quantify
#' @param pixel_header Column in pData() holding the pixel-type labels (default "sample_type")
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @return MSIobject with intensity values replaced with concentration values (ng/pixel)
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summarise_cal_levels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- create_cal_curve(cal, cal_type = "cal")
#' cal <- int2conc(cal, val_slot = "intensity", pixels = "Tissue")
#'
#' @aliases int2conc
#' @export
setMethod("int2conc", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "response", pixel_header = "sample_type", pixels = "Tissue"){

            cal_list = MSIobject@calibrationInfo@cal_list
            pixel_inds = which(pData(MSIobject)[[pixel_header]] %in% pixels)
            MSIobject = MSIobject[, pixel_inds]

            spectra(MSIobject, "conc - pg/pixel") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))
            spectra(MSIobject, "conc - pg/mm2") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            no_cal_indices = c()

            for(i in seq_len(nrow(fData(MSIobject)))){

              # Look up the calibration model by feature name (cal_list is keyed
              # by name in create_cal_curve); fall back to positional for
              # back-compatibility with older objects.
              feat_nm = fData(MSIobject)$name[i]
              eqn = if(!is.null(names(cal_list)) && feat_nm %in% names(cal_list))
                      cal_list[[feat_nm]] else cal_list[[i]]

              if(typeof(eqn) != "list"){

                message(sprintf("No calibration curve for feature %s", fData(MSIobject)$name[i]))
                no_cal_indices = c(no_cal_indices, i)

              } else{

                pixel_ints = spectraData(MSIobject)[[val_slot]][i,]
                pixel_concs = vapply(pixel_ints,
                                     function(y) chemCal::inverse.predict(eqn, y)$Prediction,
                                     numeric(1))

                mm2_scalar = (1000 / as.numeric(experimentData(MSIobject)$pixelSize)) ^2
                mm2_conc = pixel_concs * mm2_scalar

                spectra(MSIobject, "conc - pg/pixel")[i,] = pixel_concs
                spectra(MSIobject, "conc - pg/mm2")[i,] = mm2_conc

              }
            }

            if(length(no_cal_indices) > 0){
              MSIobject = MSIobject[-no_cal_indices, ]
            }

            return(MSIobject)
          })
