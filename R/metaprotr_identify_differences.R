#' identify_differences
#'
#' Shows the most differential taxonomic elements between two 
#' conditions or samples from a list defined as "spectral_count_object" 
#' with taxonomic classification. These elements are those with an 
#' absolute log2 (condition + 1 / reference + 1) > 3. If a given condition has
#' several replicates the mean value is taken into account.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts
#'    and organized by taxonomic levels. The format of this object is 
#'    similar to that generated from the function "crumble_taxonomy".
#'
#' @param target_variable Character indicating the variable name 
#'    containing the conditions or samples to be compared. This value 
#'    corresponds to the name of one column from metadata.
#'
#' @param list_conditions Atomic vector indicating two conditions to 
#'    be compared. The first element will be considered as the reference.
#'
#' @param filter_ratio Numeric value indicating the fold change filter 
#'    to be considered for the pairwise comparison. The minimal value can be 
#'    a fold change of 1.25. The default value is set at 3.
#'
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create an object in the workstation of the user.
#'
#' @return Barplots (pdf) and a csv file with the defferential taxonomic 
#'    elements between TWO conditions or sample. These elements are those that 
#'    fulfill the ratio log2 (condition + 1 / reference + 1) > filter_ratio.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(species_fw)
#' identify_differences(species_fw, "Methods", c("S", "S_EF"))
#' 
#' identify_differences(species_fw, "Methods", c("EF", "S_EF"), filter_ratio = 1.3)
#' 
#' \dontshow{setwd(.old_wd)}
#' @export

# This file is part of metaprotr.
#
# metaprotr is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free 
# Software Foundation, either version 3 of the License, or (at your option) 
# any later version. 
# metaprotr is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details. 
# You should have received a copy of the GNU General Public License along 
# with metapror. If not, see http://www.gnu.org/licenses/

identify_differences <-
  
function(spectral_count_object, target_variable, 
         list_conditions, filter_ratio = 3, force = FALSE){
  
  condition <- log2_ratio <- NULL
  spectral_count_object <- spectral_count_object
  if (!force && interactive()) {
    response <-
      select.list(
        c("yes", "no"), 
        title = "Do you allow to create a pdf and csv file containing the differential elements?")
    
    if (response == "yes"){
      if (length(spectral_count_object) == 4){
        if (spectral_count_object[[4]] == "spectral_count_object"){
          
          target_variable <- target_variable
          metadata <- spectral_count_object[[2]]
          
          if (target_variable %in% colnames(metadata) == TRUE) {
            if (all((unique(list_conditions)) %in% 
                    metadata[[target_variable]]) == TRUE) {
              if ((length(list_conditions) == 2) &
                  (length(unique(list_conditions)) == 2)) {
                
                pepsprots <- spectral_count_object[[3]]
                
                if("species.genus.family.order.class.phylum.superkingdom" %in% 
                   colnames(pepsprots)){
                  
                  if("species" %in% colnames(pepsprots) | 
                     "genus" %in% colnames(pepsprots) |
                     "family" %in% colnames(pepsprots) | 
                     "order" %in% colnames(pepsprots) | 
                     "class" %in% colnames(pepsprots) | 
                     "phylum" %in% colnames(pepsprots) | 
                     "superkingdom" %in% colnames(pepsprots)){
                    
                    if(is.numeric(filter_ratio) == TRUE &
                       abs(filter_ratio) >= 1.25 &
                       abs(filter_ratio) != Inf){
                      
                      sc_data <- spectral_count_object[[1]]
                      colnames(sc_data) <- metadata[[target_variable]]
                      sc_data <- as.data.frame(
                        sapply(split.default(sc_data, names(sc_data)), rowMeans)
                      )
                      ref <- list_conditions[1]
                      cond <- list_conditions[2]
                      ref_log <- paste(list_conditions[1], "log", sep = "_")
                      cond_log <- paste(list_conditions[2], "log", sep = "_")
                      # Use of OR instead, depending on pipeline 
                      # sc_data_pairwise <- sc_data[-which(sc_data[[ref]]==0 |
                      #   sc_data[[cond]]==0), c(ref, cond)]
                      # ref_log <- log2(sc_data_pairwise[,1])
                      # cond_log <- log2(sc_data_pairwise[,2])
                      sc_data_pairwise <- sc_data[which(sc_data[[ref]] != 0 &
                                                         sc_data[[cond]] != 0), 
                                                  c(ref, cond)]
                      sc_data_pairwise[[ref_log]] <- log2(sc_data_pairwise[, 1] + 1)
                      sc_data_pairwise[[cond_log]] <- log2(sc_data_pairwise[, 2] + 1)
                      sc_data_pairwise$log2_ratio <- sc_data_pairwise[[cond_log]] - sc_data_pairwise[[ref_log]]
                      to_plot <- sc_data_pairwise[
                        which(abs(sc_data_pairwise[, 5]) >= log2(filter_ratio)), ]
                      to_plot$condition <- ifelse(
                        to_plot$log2_ratio >= log2(filter_ratio), cond, ref
                      )
                      if(dim(to_plot)[1] == 0){
                        stop(paste("No taxonomic entities were found above a ratio of", filter_ratio, sep = " "))
                      }
                      else{
                        # taxonomic level
                        if("species" %in% colnames(pepsprots)){
                          taxa <- "Species"
                        }
                        else if("genus" %in% colnames(pepsprots)){
                          taxa <- "Genera"
                        }
                        else if("family" %in% colnames(pepsprots)){
                          taxa <- "Families"
                        }
                        else if("order" %in% colnames(pepsprots)){
                          taxa <- "Orders"
                        }
                        else if("class" %in% colnames(pepsprots)){
                          taxa <- "Classes"
                        }
                        else if("phylum" %in% colnames(pepsprots)){
                          taxa <- "Phylum"
                        }
                        else{
                          taxa <- "Superkingdoms"
                        }
                        
                        title <- paste(taxa, " abundances ",
                                       list_conditions[2], "/", list_conditions[1], 
                                       " (spectral counts ratio)", sep = ""
                        )
                        xlab <- paste("Log2(", cond, "/", ref, ")", sep = " ")
                        to_plot[[taxa]] <- row.names(to_plot)
                        # CSV file
                        filename <- paste("differential", taxa, cond, "vs", ref, 
                                          filter_ratio, sep = "_")
                        f_csv <- paste(filename, ".csv", sep = "")
                        write.table(to_plot, f_csv, append = FALSE, sep = "\t", 
                                    dec = ".", row.names = FALSE, col.names = TRUE)
                        print(paste(f_csv, "was generated", sep = " "))
                        
                        # Plot
                        # To keep graphic settings
                        old_par <- par(no.readonly = TRUE)
                        on.exit(suppressWarnings(par(old_par)))
                        
                        to_plot$taxa <- to_plot[[taxa]]
                        p <- ggplot(to_plot, aes(x = reorder(taxa, log2_ratio), 
                                                 y = log2_ratio, 
                                                 color = condition)) +
                          geom_bar(stat = "identity", fill = "white") +
                          coord_flip() +
                          labs(title = title, y = xlab, x = " ") +
                          theme_minimal(base_size = 11)
                        f_pdf <- paste(filename, ".pdf", sep = "")
                        # for aestethics
                        height_plot <- dim(to_plot)[1] * 0.15
                        if(height_plot < 4){
                          height_plot <- 4
                        }
                        else if(height_plot > 21){
                          height_plot <- 19
                        }
                        else{
                        }
                        pdf(f_pdf, width = 9, height = height_plot)
                        print(p)
                        dev.off()
                        print(paste(f_pdf, "was generated", sep = " "))
                      }
                    }
                    else{
                      stop("The third argument must be higher than 1.25")
                    }
                  }
                  else{
                    stop("Spectral counts must be reorganized with the function crumble_taxonomy")
                  }
                }
                else{
                  stop("Taxonomy must be added with the function add_taxonomy")
                }
              }
              else{
                stop("Precise only TWO conditions in the third argument")
              }
            }
            
            else{
              vars <- paste(list_conditions, collapse = ", ")
              stop(paste(
                c("The variables [", vars, "] must be present in: '", 
                  target_variable, "'"), 
                collapse = " ")
              )
            }
          }
          
          else{
            options <- colnames(metadata)
            stop(paste(
              c("Change target_variable for ONE of these options: ", options), 
              collapse = "' '")
            )
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
