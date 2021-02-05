#' plot_intensities
#'
#' Draws violin plots containing the abundance intensities expressed as 
#' spectral counts per level (peptides, subgroups, groups or taxonomic entities) 
#' in provided samples or conditions from a "spectral_count_object". If the 
#' provided conditions have several replicates the mean value is taken into account.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts 
#'    from peptides, subgroups, groups or taxonomic levels. The format of 
#'    this object is similar to that generated from the functions 
#'    "getsc_specific" and "crumble_taxonomy".
#'    
#' @param target_variable Character indicating the name of one column from 
#'    metadata, the column must contain the conditions to be displayed.
#' 
#' @param image_title Character indicating the title to be displayed 
#'    in the generated image.
#'
#' @param force Logic value set at FALSE by default in order to ask
#'    for permission to create a pdf file in the workstation of the user. 
#'    
#' @return Violin plots (pdf) indicating the spectral counts of the different 
#'    levels (peptides, subgroups, groups or taxonomic entities) per sample 
#'    or condition.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' plot_intensities(fecal_waters, "SC_name", "Title to display inside the plot")
#' 
#' data(species_fw)
#' plot_intensities(species_fw, "Condition", "Abundance per condition")
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

plot_intensities <-
  
function(spectral_count_object, target_variable, image_title = NULL, force = FALSE){
  
  spectra <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <- select.list(
      c("yes", "no"),
      title = paste("Do you allow to create a violin plot with all abundances per", target_variable, "?", sep = " ")
    )
    
    if (response == "yes") {
      if (length(spectral_count_object) == 4) {
        if (spectral_count_object[[4]] == "spectral_count_object") {
          
          target_variable <- target_variable
          metadata <- spectral_count_object[[2]]
          
          if (target_variable %in% colnames(metadata) == TRUE) {
            if(length(unique(metadata[[target_variable]])) <= 1){
              stop("Check that there are at least two conditions in the provided data")
            }
            else{
              if (is.character(image_title) == TRUE && nchar(image_title) > 1){
                
                sc_data <- spectral_count_object[[1]]
                colnames(sc_data) <- metadata[[target_variable]]
                # aggregate columns by the mean of the same condition
                sc_data <- sapply(split.default(sc_data, names(sc_data)), rowMeans)
                sc_data <- as.data.frame(sc_data)
                sc_data_stacked <- stack(sc_data[, 1:length(colnames(sc_data))])
                names(sc_data_stacked)[1] <- "spectra"
                names(sc_data_stacked)[2] <- target_variable
                sc_data_stacked <- cbind.data.frame(sc_data_stacked, element = rep(rownames(sc_data), length(colnames(sc_data))))
                sc_data_stacked$spectra <- log2(sc_data_stacked$spectra + 1)
                
                # Plot
                if (names(spectral_count_object[1]) == "SC_subgroups") {
                  elements <- "Metaprotein"
                }
                else if (names(spectral_count_object[1]) == "SC_groups") {
                  elements <- "Metagroup"
                }
                else{
                  elements <- "Peptide"
                }
                main <- image_title
                xlab <- ""
                ylab <- "log2 (spectral count + 1)"
                filename <- paste(elements, "intensities_by", target_variable, length(colnames(sc_data)), "elements", sep = "_")
                filename <- paste(filename, ".pdf", sep = "")
                variable_x <- sc_data_stacked[[target_variable]]
                
                # To keep graphics settings
                old_par <- par(no.readonly = TRUE)
                on.exit(suppressWarnings(par(old_par)))
                
                p <- ggplot(
                  sc_data_stacked, aes(x = variable_x, y = spectra, fill = variable_x)) +
                  geom_violin(trim = FALSE) +
                  geom_boxplot(width = 0.03, fill = "white") +
                  scale_fill_viridis_d(direction = -1) +
                  labs(title = main, x = xlab, y = ylab) +
                  theme_classic(base_size = 19) +
                  theme(legend.position = "none")
                if( length(colnames(sc_data)) > 15 ){
                  width <- length(colnames(sc_data)) * 0.35
                }
                else{
                  width <- 12
                }
                pdf(filename, width = width, height = 10)
                print(p)
                dev.off()
                print(paste("Graph", filename, "was created!", sep = " "))
              }
              else{
                stop("The third argument must be a character indicating the name to diplay inside the plot")
              }
            }
          }
          else{
            options <- colnames(metadata)
            stop(paste(c("Change target_variable for ONE of these options: ", options), collapse = "' '"))
          }
        }
        else {
          stop("Invalid spectral count object")
        }
      }
      else{
        stop("Invalid spectral count object")
      }
    }
    
    else{
      stop("No file was created")
    }
  }
}
