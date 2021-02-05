#' inspect_sample_elements
#'
#' Displays a graph that indicates the number of common elements from a 
#' "spectral_count_object" (peptides, subgroups, groups or taxonomic entities) 
#' per sample. This function is useful to distinguish heterogeneity between 
#' samples in an experimental design.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts. 
#'    The spectral data can be organized by peptides, subgroups, 
#'    groups or taxonomic levels.
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a pdf file in the workstation of the user. 
#' 
#' @return Barplots (pdf) ilustrating the common spectral elements 
#'    (peptides, subgroups, groups, taxonomic elements) per sample
#'    in a "spectral_count_object".
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' inspect_sample_elements(fecal_waters)
#' 
#' \dontshow{setwd(.old_wd)}
#' @export

# This file is part of metaprotr.
#
# metaprotr is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation, 
# either version 3 of the License, or (at your option) any later version. 
# metaprotr is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details. 
# You should have received a copy of the GNU General Public License along with metapror. 
# If not, see http://www.gnu.org/licenses/


inspect_sample_elements <- 

function(spectral_count_object, force = FALSE){
  
  counter <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <-
      select.list(c("yes", "no"), 
                  title = "Allow to create a pdf file with frequency of elements per sample?")
    if (response == "yes"){
      if(length(spectral_count_object) == 4){
        if(spectral_count_object[[4]] == "spectral_count_object"){
          
          #---------  Check type of spectral_count_object for plot info
          peps_prots <- spectral_count_object[[3]]
          if("species" %in% colnames(peps_prots) | "genus" %in% colnames(peps_prots) | 
             "family" %in% colnames(peps_prots) | "order" %in% colnames(peps_prots) | 
             "class" %in% colnames(peps_prots) | "phylum" %in% colnames(peps_prots) | 
             "superkingdom" %in% colnames(peps_prots)){
            spectral_origin <- colnames(peps_prots)[length(colnames(peps_prots))]
          }
          else{
            type_sc_data <- names(spectral_count_object[1])
            if(type_sc_data == "SC_specific_peptides"){
              spectral_origin <- "peptides"
            }
            else if(type_sc_data == "SC_subgroups"){
              spectral_origin <- "subgroups"
            }
            else{
              spectral_origin <- "groups"
            }
          }
          #---------  Count elements' presence across the different samples
          sc_data <- spectral_count_object[[1]]
          decision_matrix <- stack(sc_data[, 1:length(colnames(sc_data))])
          names(decision_matrix)[1] <- "spectra"
          names(decision_matrix)[2] <- "sample"
          decision_matrix <- cbind.data.frame(
            decision_matrix, 
            element = rep(rownames(sc_data), length(colnames(sc_data))))
          decision_matrix$presence <- ifelse(decision_matrix$spectra == 0, 0, 1) 
          decision_matrix <- tapply(decision_matrix$presence, 
                                    list(decision_matrix$element, decision_matrix$sample), 
                                    FUN=sum)
          decision_matrix <- as.data.frame(decision_matrix)
          elements_count <- decision_matrix[, 1:length(colnames(sc_data))] > 0 
          elements_count <- apply(elements_count, 1, sum)
          decision_matrix$elements <- as.factor(row.names(decision_matrix))
          decision_matrix$counter <- elements_count
          
          #---------  Output
          
          # To keep graphic settings
          old_par <- par(no.readonly = TRUE)
          on.exit(suppressWarnings(par(old_par)))
          
          ylab = paste("Number of", spectral_origin, sep = " ")
          p <- ggplot(decision_matrix, aes(x = counter)) + 
               geom_histogram(fill = "darkorange", bins = length(colnames(sc_data)) * 2) +
               scale_fill_viridis_d(direction = -1) +
               scale_x_continuous(breaks=c(1:length(colnames(sc_data)))) +
               labs(title = "Elements shared by N samples", x = "Samples", y = ylab) +
               theme_classic(base_size = 15) 
          filename <- paste("frequency", spectral_origin, "per_sample", sep = "_")
          filename <- paste(filename, ".pdf", sep = "")
          if(length(colnames(sc_data)) > 15 ){
            width = length(colnames(sc_data)) * 0.35
          }
          else{
            width = 11
          }
          pdf(filename, width = width, height = 9)
          print(p)
          dev.off()
          print(paste("Graph", filename, "was created!", sep = " "))
        }
        else{
          stop("Invalid object")
        }
      }
      else{
        stop("Invalid object")
      }
    }
    else{
      stop("No file was created")
    }
  }
}