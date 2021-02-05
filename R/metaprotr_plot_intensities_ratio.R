#' plot_intensities_ratio
#'
#' Generates a scatter plot of the log2 (ratio + 1) between two conditions considering
#' the spectral counts of each entity (peptides, subgroups, groups or taxonomic 
#' levels) from a "spectral_count_object". If a given condition has several 
#' replicates the mean value is taken into account.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts
#'    organized by peptides, subgroups, groups or taxonomic levels. 
#'    The format of this object is similar to that generated from 
#'    the functions "getsc_specific" and "crumble_taxonomy".
#'
#' @param target_variable Character indicating the variable name 
#'    containing the conditions to be compared. This value corresponds 
#'    to the name of one column from metadata.
#'
#' @param list_conditions Atomic vector indicating two conditions to A
#'    be compared. The first element is considered as the reference (denominator) 
#'    for ratio calculations.
#'
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create object in the workstation of the user.
#'
#' @return A scatter plot (pdf) indicating the log2 (ratio + 1) of the entities
#'    (peptides, soubgroups, groups, taxonomic) between the two conditions provided.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' 
#' plot_intensities_ratio(fecal_waters, "Methods", c("EF", "S"))
#' 
#' plot_intensities_ratio(fecal_waters, "SC_name", c("Q1", "Q2"))
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

plot_intensities_ratio <-
  
function(spectral_count_object, target_variable, list_conditions, force = FALSE){
  
  spectral_count_object <- spectral_count_object
  if (!force && interactive()) {
    response <-
      select.list(c("yes", "no"), title = "Do you allow to create a pdf file with log2 abundance ratios?")
    
    if (response == "yes") {
      if (length(spectral_count_object) == 4) {
        if (spectral_count_object[[4]] == "spectral_count_object") {
          
          target_variable <- target_variable
          metadata <- spectral_count_object[[2]]
          
          if (target_variable %in% colnames(metadata) == TRUE) {
            if (all((unique(list_conditions)) %in% metadata[[target_variable]]) == TRUE) {
              if ((length(list_conditions) == 2) &
                  (length(unique(list_conditions)) == 2)) {
                
                sc_data <- spectral_count_object[[1]]
                colnames(sc_data) <- metadata[[target_variable]]
                
                # aggregate data in case the used selected a condition comparison that has replicates
                sc_data <- as.data.frame(sapply(split.default(sc_data, names(sc_data)), rowMeans))
                ref <- list_conditions[1]
                cond <- list_conditions[2]
                
                # Use of OR instead, depending of pipeline 
                # sc_data_pairwise <- sc_data[-which(sc_data[[ref]]==0 | sc_data[[cond]]==0), c(ref, cond)]
                # ref_log <- log2(sc_data_pairwise[,1])
                # cond_log <- log2(sc_data_pairwise[,2])
                sc_data_pairwise <- sc_data[-which(sc_data[[ref]] == 0 & sc_data[[cond]] == 0), c(ref, cond)]
                ref_log <- log2(sc_data_pairwise[, 1] + 1)
                cond_log <- log2(sc_data_pairwise[, 2] + 1)
                
                # Matrix creation containing the abundances in log2 and ratio = cond/ref in descendent order (by reference)
                to_plot <- matrix(NA, length(ref_log), 3)
                to_plot[, 1] <- ref_log
                to_plot[, 2] <- cond_log
                to_plot[, 3] <- cond_log - ref_log
                to_plot <- to_plot[rev(order(to_plot[, 1])),]
                
                # Plot logic with marks at |log2(ratio)| > 3
                pepsprots <- spectral_count_object[[3]]
                # For labels in plot
                if("species" %in% colnames(pepsprots)){
                  taxa <- "Species"
                }
                else if("genus" %in% colnames(pepsprots)){
                  taxa <- "Genus"
                }
                else if("family" %in% colnames(pepsprots)){
                  taxa <- "Family"
                }
                else if("order" %in% colnames(pepsprots)){
                  taxa <- "Order"
                }
                else if("class" %in% colnames(pepsprots)){
                  taxa <- "Class"
                }
                else if("phylum" %in% colnames(pepsprots)){
                  taxa <- "Phylum"
                }
                else if("superkingdom" %in% colnames(pepsprots)){
                  taxa <- "Superkingdom"
                }
                else{
                  taxa <- ""
                }
                if (names(spectral_count_object[1]) == "SC_subgroups") {
                  if (nchar(taxa) == 0){
                    elements <- "Subgroup"
                  }
                  else{
                    elements <- paste(c(taxa, "(from subgroups)"), collapse = " ")
                  }
                }
                else if (names(spectral_count_object[1]) == "SC_groups") {
                  if (nchar(taxa) == 0){
                    elements <- "Group"
                  }
                  else{
                    elements <- paste(c(taxa, "(from groups)"), collapse = " ")
                  }
                }
                else{
                  if (nchar(taxa) == 0){
                    elements <- "Peptide"
                  }
                  else{
                    elements <- paste(c(taxa, "(from peptides)"), collapse = " ")
                  }
                }
                number <- paste(c("n=", length(ref_log)), collapse = "")
                xlab <- paste(c(elements, " by decreasing abundance in ", ref, " (", number, ")"), collapse = "")
                ylab <- paste(c("log2(ratio ", cond, "/", ref, " )"), collapse = "")
                main <- paste(c("Abudances in ", cond, " versus ", ref), collapse = "")
                filename <- paste(elements, "ratios", cond, "VS", ref, sep = "_")
                filename <- paste(filename, ".pdf", sep = "")
                
                # To keep graphics settings
                old_par <- par(no.readonly = TRUE)
                on.exit(suppressWarnings(par(old_par)))
                
                pdf(filename, width = 7.5, height = 6)
                par(fig = c(0, 0.85, 0, 1))
                plot(to_plot[, 3], pch = ".", cex = 2, xlab = xlab, ylab = ylab, main = main)
                abline(h = 0, col = "blue", lwd = 0.5)
                abline(h = -1.584963, col = "blue", lwd = 1)
                abline(h = 1.584963, col = "blue", lwd = 1)
                
                par(fig = c(0.70, 1, 0, 1), new = TRUE)
                boxplot(to_plot[, 3], axes = FALSE)
                abline(h = 0, col = "blue", lwd = 0.5)
                abline(h = -1.584963, col = "blue", lwd = 1)
                abline(h = 1.584963, col = "blue", lwd = 1)
                dev.off()
                print(paste("The figure", filename, "was generated", sep = " "))
              }
              
              else{
                stop("Precise only TWO conditions in list_conditions argument")
              }
            }
            
            else{
              vars <- paste(list_conditions, collapse = ", ")
              stop(paste(c("The variables [", vars, "] must be present in: '", target_variable, "'"), collapse = " "))
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
