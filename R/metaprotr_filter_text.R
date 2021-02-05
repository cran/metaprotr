#' filter_text
#'
#' Matches the entities containing a given chain of characters inside an explanatory 
#' variable (column name) of the dataframe "peptides_proteins" from a 
#' "spectral_count_object". Based on the user's decision, the peptides, subgroups, 
#' groups or taxonomic levels containig the provided chain of characters will 
#' be kept or discarted in a newly generated object. 
#' 
#' @param spectral_count_object List containing dataframes with proteomics
#'    elements whose abundance is expressed as spectral counts and are
#'    organized by peptides, subgroups, groups or taxonomic levels. 
#'    The format of this object is similar to that generated from 
#'    the functions "getsc_specific" and "crumble_taxonomy".
#' 
#' @param pepsprots_feature Character indicating the name of one explanatory 
#'    variable (ONE column name) of the dataframe "peptides_proteins". 
#' 
#' @param text_to_filter Character containig the text to be searched in 
#'    the "pepsprots_feature" content.
#' 
#' @param decision Character indicating wether the elements containing 
#'    the matched text will be kept or dirscarted. The two allowed option 
#'    are: "keep" or "discard".
#' 
#' @return A list defined as "spectral_count_object" with or without the 
#'     elements (peptides, subgroups, groups, taxonomic items) that 
#'     matched the provided text in a given variable of the 
#'     "peptides_proteins" dataframe.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' data(species_fw)
#' 
#' cysteine_alkylations <- filter_text(fecal_waters, "Modifs", "57.02146", "keep")
#' 
#' exclude_medimonas <- filter_text(species_fw, "organism", "Merdimonas faecis BR31", "discard")
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


filter_text <- 

function(spectral_count_object, pepsprots_feature, text_to_filter, decision){
  
  spectral_count_object <- spectral_count_object
  
  if(length(spectral_count_object) == 4){
    if(spectral_count_object[[4]] == "spectral_count_object"){
      
      type_sc_data <- names(spectral_count_object[1])
      pepsprots_feature <- pepsprots_feature
      metadata <- spectral_count_object[[2]]
      peps_prots <- spectral_count_object[[3]]

      if(is.character(pepsprots_feature) == TRUE && 
         pepsprots_feature %in% colnames(peps_prots) == TRUE){
        
        peps_prots[[pepsprots_feature]] <- as.character(peps_prots[[pepsprots_feature]])
        
        if(is.character(text_to_filter) == TRUE && 
           nchar(text_to_filter) > 0){
          
          if(decision == "keep" | decision == "discard"){
            #------------- Match text in the provided feature 
            if(type_sc_data == "SC_specific_peptides"){
              spectral_origin <- "Peptide"
              selected_elements <- peps_prots[
                grep(text_to_filter, peps_prots[[pepsprots_feature]]), "Peptide"]
            }
            else if(type_sc_data == "SC_subgroups"){
              spectral_origin <- "Subgroup"
              selected_elements <- peps_prots[
                grep(text_to_filter, peps_prots[[pepsprots_feature]]), "Subgroup"]
            }
            else if (type_sc_data == "SC_groups"){
              spectral_origin <- "Group"
              selected_elements <- peps_prots[
                grep(text_to_filter, peps_prots[[pepsprots_feature]]), "Group"]
            }
            else{
              stop("Invalid spectral_count_object")
            }
            
            #------------- Elements to keep or discard
            if(decision == "keep"){
              new_peps_prots <- peps_prots[which(peps_prots[[spectral_origin]] %in% selected_elements), ]
            }
            else{
              new_peps_prots <- peps_prots[-which(peps_prots[[spectral_origin]] %in% selected_elements), ]
            }
            # To point at spectral data formated by taxonomy
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
            }
            elements_to_keep <- as.vector(new_peps_prots[[spectral_origin]])
            sc_data <- spectral_count_object[[1]]
            new_sc_data <- sc_data[rownames(sc_data) %in% elements_to_keep, ]
            
            #------------- Status message
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
            stop("The fourth argument must be either 'keep' or 'discard'")
          }
        }
        else{
          stop("The third argument must be a chain of characters")
        }
      }
      else{
        options_pepsprots <- paste(colnames(peps_prots), collapse = "', '")
        stop(paste(c("The second argument is invalid. This argument must be equal to one these options: '", options_pepsprots), collapse = " "))
      }
    }
    else{
      stop("Invalid spectral count object")
    }
  }
  else{
    stop("Invalid spectral count object")
  }
}