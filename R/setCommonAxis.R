setGeneric("setCommonAxis", function(MSIobjects, ...) standardGeneric("setCommonAxis"))

#' Align a collection of MSI objects to a common feature axis
#'
#' Puts every object in a list onto the same feature set, in the same order, so
#' they can be combined. Features missing from an object are padded, which is
#' what makes acquisitions run on slightly different MRM panels comparable.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobjects A list of `MSImagingExperiment` objects.
#' @param ref_fdata A `MassDataFrame` giving the reference feature axis, e.g.
#'   `fData()` of whichever object carries the full panel.
#' @return A list of the same objects, each carrying `ref_fdata` as its
#'   `fData()` with spectra reordered and padded to match.
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' objs <- list(readRDS(p1), readRDS(p2))
#' aligned <- setCommonAxis(objs, fData(objs[[1]]))
#'
#' @family combining acquisitions
#' @aliases setCommonAxis
#' @export
setMethod("setCommonAxis", "list",
          function(MSIobjects, ref_fdata){

            ref_features = data.frame(ref_fdata)

            # Common axis
            for(i in seq_along(MSIobjects)){

              MSIobject = MSIobjects[[i]]

              # if no ref feature in MSIobject
              if(! any( ref_features$name %in% fData(MSIobject)$name) ){
                MSIobjects[[i]] = NULL
              }


              # Remove those features not in reference
              MSIobject = MSIobject[which(fData(MSIobject)$name %in% ref_features$name), ]



              # Add blank channels for ref features not in MSIobject
              blank_feat_add = ref_features$name[which(! ref_features$name %in% fData(MSIobject)$name)]
              new_fdata = data.frame(fData(MSIobject))
              new_fdata$mz = seq_len(nrow(new_fdata))

              new_idata = spectra(MSIobject)

              for(feat_name in blank_feat_add){
                # Index of new position in MSIobject
                ind = nrow(new_fdata) +1

                # Update feature metadata
                new_line = ref_features[which(ref_features$name == feat_name), ]
                new_line$mz = ind
                new_fdata = rbind(new_fdata, new_line)

                # Add blank channel into spectra
                new_idata = rbind(new_idata, matrix(nrow=1, ncol=ncol(MSIobject)))
              }

              rownames(new_idata) = new_fdata$mz

              # Correct the order
              out_idata = matrix(NA, nrow=nrow(ref_features), ncol = ncol(MSIobject))
              for(row_ind in seq_len(nrow(ref_features))){

                feat = ref_features$name[row_ind]

                original_ind = which(new_fdata$name == feat)
                out_idata[row_ind, ] = new_idata[original_ind, ]

              }

              MSIobjects[[i]] <- MSImagingExperiment(spectraData=out_idata,
                                                     featureData=ref_fdata,
                                                     pixelData=pData(MSIobject),
                                                     experimentData = experimentData(MSIobject))
            }

            return(MSIobjects)

          })
