#' filter_unshared
#'
#' Keeps the elements from a "spectral_count_object" that are specific to one
#' level of a provided explanatory variable. This variable MUST correspond to
#' the name of ONE column of the dataframe "metadata" (eg. conditions or samples).
#' This function allows to identify the elements that are specific to one 
#' condition or one sample. 
#' 
#' @param spectral_count_object List containing dataframes with proteomics
#'    elements whose abundance is expressed as spectral counts and are
#'    organized by peptides, subgroups, groups or taxonomic levels. 
#'    The format of this object is similar to that generated from 
#'    the functions "getsc_specific" and "crumble_taxonomy".
#' 
#' @param metadata_feature Character indicating the name of one explanatory 
#'    variable (ONE column name) of the dataframe "metadata".
#' 
#' @return A list defined as "spectral_count_object" with the specific 
#'    elements per sample or condition, having at least one spectra.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' data(species_fw)
#' 
#' specific_elements_per_sample <- filter_unshared(fecal_waters, "SampleID")
#' 
#' specific_elements_per_condition <- filter_unshared(species_fw, "Condition")
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

filter_unshared <- 

function(spectral_count_object, metadata_feature){
  
  spectral_count_object <- spectral_count_object
  
  if(length(spectral_count_object) == 4){
    if(spectral_count_object[[4]] == "spectral_count_object"){
      
      type_sc_data <- names(spectral_count_object[1])
      metadata <- spectral_count_object[[2]]
      
      if(metadata_feature %in% colnames(metadata) == TRUE){
        
        #---------  Agregate Spectral counts in function of the provided feature 
        sc_data <- spectral_count_object[[1]]
        colnames(sc_data) <- metadata[[metadata_feature]]
        sc_data <- sapply(split.default(sc_data, names(sc_data)), rowSums)
        sc_data <- as.data.frame(sc_data)
        
        #---------  Check that the feature has specific elements
        decision_matrix <- stack(sc_data[, 1:length(colnames(sc_data))])
        names(decision_matrix)[1] <- "spectra"
        names(decision_matrix)[2] <- "feature"
        decision_matrix <- cbind.data.frame(
          decision_matrix, 
          element = rep(rownames(sc_data), length(colnames(sc_data)))
        )
        decision_matrix$presence <- ifelse(decision_matrix$spectra == 0, 0, 1) 
        decision_matrix <- tapply(decision_matrix$presence, 
                                  list(decision_matrix$element, decision_matrix$feature), 
                                  FUN=sum)
        decision_matrix <- as.data.frame(decision_matrix)
        elements_count <- decision_matrix[, 1:length(colnames(sc_data))] > 0 
        elements_count <- apply(elements_count, 1, sum)
        decision_matrix$elements <- as.factor(row.names(decision_matrix))
        decision_matrix$counter <- elements_count
        unique_elements <- as.vector(decision_matrix$elements[decision_matrix$counter == 1])
        new_sc_data <- spectral_count_object[[1]]
        new_sc_data <- new_sc_data[rownames(new_sc_data) %in% unique_elements, ]
        
        #---------  Update rows in peps_prots
        peps_prots <- spectral_count_object[[3]]
        if("species" %in% colnames(peps_prots) | "genus" %in% colnames(peps_prots) | 
           "family" %in% colnames(peps_prots) | "order" %in% colnames(peps_prots) | 
           "class" %in% colnames(peps_prots) | "phylum" %in% colnames(peps_prots) | 
           "superkingdom" %in% colnames(peps_prots)){
          if("species" %in% colnames(peps_prots)){
            spectral_origin <- "species"
          }
          else if("genus" %in% colnames(peps_prots)){
            spectral_origin <- "genus"
          }
          else if("family" %in% colnames(peps_prots)){
            spectral_origin <- "family"
            }
          else if("order" %in% colnames(peps_prots)){
            spectral_origin <- "order"
          }
          else if("class" %in% colnames(peps_prots)){
            spectral_origin <- "class"
          }
          else if("phylum" %in% colnames(peps_prots)){
            spectral_origin <- "phylum"
          }
          else{
            spectral_origin <- "superkingdom"
          }
          new_peps_prots <- peps_prots[peps_prots[[spectral_origin]] %in% unique_elements, ]
        }
        else{
          if(type_sc_data == "SC_specific_peptides"){
            new_peps_prots <- peps_prots[peps_prots$Peptide %in% unique_elements, ]
          }
          else if(type_sc_data == "SC_subgroups"){
            new_peps_prots <- peps_prots[peps_prots$Subgroup %in% unique_elements, ]
          }
          else{
            new_peps_prots <- peps_prots[peps_prots$Group %in% unique_elements, ]
          }
        }
        #---------  Output: status message
        print("========= Object information ============")
        print(paste("Number of samples: ", length(unique(metadata$msrunfile)), sep = ""))
        print(paste("Number of groups: ", length(unique(new_peps_prots$Group)), sep = ""))
        print(paste("Number of subgroups or metaproteins: ", length(unique(new_peps_prots$Subgroup)), sep = ""))
        print(paste("Number of peptides: ", length(unique(new_peps_prots$Peptide)), sep = ""))
        print("=========================================")
        output <- list(new_sc_data, metadata, new_peps_prots, "spectral_count_object")
        names(output)[1] <- type_sc_data
        names(output)[2] <- "metadata"
        names(output)[3] <- "peptides_proteins"
        names(output)[4] <- "type_object"
        print("Spectral count object generated")
        return(output)
      }
      else{
        options_meta <- paste(colnames(metadata), collapse = "', '")
        stop(paste(c("The second argument is invalid. This argument must be one these options: '", options_meta), collapse = " "))
      }
    }
    else{
      stop("Invalid object")
    }
  }
  else{
    stop("Invalid object")
  }
}