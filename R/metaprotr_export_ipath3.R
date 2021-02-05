#' export_ipath3
#'
#' Exports the KEGG Orthology (KO) terms in the adapted format to be used in 
#' the tool \href{https://pathways.embl.de/}{iPATH3}. The exported data is 
#' obtained from a "spectral_count_object" containing the functional 
#' annotation of the identified proteins from one condition or sample.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing protein abundance expressed as spectral counts by a 
#'    taxonomic level. The functional annotation must be added to this
#'    object. The format of this object is similar to that generated from 
#'    the function "add_kegg".
#'
#' @param type_export Character indicating the type of export to be used. The
#'    possible options are: i) "all" that selects all the KO terms from a given 
#'    sample or a given condition; and, ii) "selection" that extracts the 
#'    KO terms present in selected taxonomic entities (one or more).
#' 
#' @param target_variable Character indicating the column name from metadata 
#'    containing the condition or sample to be analyzed.
#'
#' @param sample_condition Atomic vector indicating the sample from which 
#'    the functional information will be extracted.
#' 
#' @param taxonomic_levels Optional vector indicating the taxonomic levels
#'    from which the KO terms will be extracted. This option is needed only
#'    if the type of export is "selection".
#' 
#' @param hexadecimal_color Character indicating the color to be used in 
#'    iPATH3, this value must be indicated in hexadecimal format (eg. #ff0000).
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a file in the workstation of the user.
#'
#' @return A csv file containing the KO terms present in a given sample or 
#'    condition. The content of this file can be inserted directly in the tool  
#'    \href{https://pathways.embl.de/}{iPATH3}. The width of the lines in 
#'    iPATH3 will be displayed by the  percentage of spectra in the selected 
#'    sample or condition. In this way, KO terms belonging to a given taxonomic 
#'    level are represented in three intervals based on their abundace: i) 
#'    below 2 percent, i) between 2 to 10 percent, or i) above 10 percent.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(species_annot_fw)
#' 
#' export_ipath3(
#'    species_annot_fw, 
#'    "all", 
#'    "SampleID", 
#'    "Q1_prot", 
#'    "#840AA3"
#' )
#' 
#' taxonomic_entities <- c("Bacteroides caccae", "Coprococcus catus", "Merdimonas faecis")
#' export_ipath3(
#'    species_annot_fw, 
#'    "selection", 
#'    "SC_name", 
#'    "FW2", 
#'    "#28c1df", 
#'    taxonomic_entities
#' )
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

export_ipath3 <-
  
function(spectral_count_object, type_export, target_variable, 
         sample_condition, hexadecimal_color, 
         taxonomic_levels = NULL, force = FALSE){
  
  ko <- percent_spectra <- NULL
  spectral_count_object <- spectral_count_object
  
  if(!force && interactive()){
    response <- select.list(
      c("yes", "no"), 
      title = "Do you allow to create a csv file with KO terms?"
    )
    if(response == "yes"){
      if(length(spectral_count_object) == 4){
        if(spectral_count_object[[4]] == "spectral_count_object"){
          
          pepsprots <- spectral_count_object[[3]]
          if("ko" %in% colnames(pepsprots)){
            
            metadata <- spectral_count_object[[2]]
            if(target_variable %in% colnames(metadata)){
              
              options_var <- metadata[[target_variable]]
              if(sample_condition %in% options_var){
                
                if("species" %in% colnames(pepsprots)){
                  taxa <- "species"
                }
                else if("genus" %in% colnames(pepsprots)){
                  taxa <- "genus"
                }
                else if("family" %in% colnames(pepsprots)){
                  taxa <- "family"
                }
                else if("order" %in% colnames(pepsprots)){
                  taxa <- "order"
                }
                else if("class" %in% colnames(pepsprots)){
                  taxa <- "class"
                }
                else if("phylum" %in% colnames(pepsprots)){
                  taxa <- "phylum"
                }
                else{
                  taxa <- "superkingdom"
                }
                
                sc_data <- spectral_count_object[[1]]
                colnames(sc_data) <- metadata[[target_variable]]
                sc_data <- round(sapply(split.default(sc_data, names(sc_data)), rowMeans), 2)
                sc_data <- as.data.frame(sc_data)
                sc_data <- sc_data[which(sc_data[[sample_condition]] > 0), sample_condition, drop = FALSE]
                total_spectra <- sum(sc_data[[sample_condition]])
                sc_data$percent_spectra <- 100 * (sc_data[[sample_condition]] / total_spectra)
                sc_data[[taxa]] <- rownames(sc_data)
                
                if(type_export == "all"){
                  
                  names(sc_data)[1] <- "spectra"
                  sc_data <- merge(sc_data, pepsprots[, c(taxa, "ko")], by = taxa, all.x = TRUE)
                  ipath_format <- sc_data %>% group_by(ko) %>% summarize(percent_spectra = max(percent_spectra))
                  ipath_format$width <- ifelse(
                    ipath_format$percent_spectra > 10, "W25", 
                    ifelse(ipath_format$percent_spectra < 2, "W5", "W13")
                  )
                  ipath_format$color <- hexadecimal_color
                  ipath_format <- ipath_format[, c("ko", "color", "width")]
                  filename <- paste("all", taxa, target_variable, sample_condition, "ipath3", sep = "_")
                  filename <- paste(filename, ".csv", sep = "")
                  write.table(ipath_format, filename, append = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
                  print(paste(filename, "was generated", sep = " "))
                  
                }
                else if(type_export == "selection"){
                  
                  if(is.null(taxonomic_levels)){
                    stop("The sixth argument can not be empty")
                  }
                  else{
                    if(all(unique(taxonomic_levels) %in% pepsprots[[taxa]]) == TRUE){
                      sc_data <- sc_data[which(sc_data[[taxa]] %in% taxonomic_levels),]
                      names(sc_data)[1] <- "spectra"
                      sc_data <- merge(sc_data, pepsprots[, c(taxa, "ko")], by = taxa, all.x = TRUE)
                      ipath_format <- sc_data %>% group_by(ko) %>% summarize(percent_spectra = max(percent_spectra))
                      ipath_format$width <- ifelse(
                        ipath_format$percent_spectra > 10, "W25", 
                        ifelse(ipath_format$percent_spectra < 2, "W5", "W13")
                      )
                      ipath_format$color <- hexadecimal_color
                      ipath_format <- ipath_format[, c("ko", "color", "width")]
                      filename <- paste("selection", taxa, target_variable, sample_condition, "ipath3", sep = "_")
                      filename <- paste(filename, ".csv", sep = "")
                      write.table(ipath_format, filename, append = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
                      print(paste(filename, "was generated", sep = " "))
                    }
                    else{
                      stop("Provide taxonomic levels present in the spectral count object")
                    }
                  }
                }
                else{
                  stop("The possible options for the second argument are: 'all' or 'selection'")
                }
              }
              else{
                vars <- unique(options_var)
                stop(paste(c("The fourth argument must be one these options", vars), collapse = ", "))
              }
            }
            else{
              var_options <- colnames(metadata)
              stop(paste(c("The third argument must be one of the following options", var_options), collapse = ", "))
            }
          }
          else{
            stop("Functional annotation must be added previously with the function add_kegg")
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