#' remove_element
#'
#' Removes elements from a "spectral_count_object". These elements can 
#' be: i) samples, ii) peptides, iii) proteins, iv) soubgroups, v) groups, 
#' vi) sequences, vii) species, viii) genus, ix) family, x) order, 
#' xi) class, xii) phylum or xiii) superkingdom. 
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts 
#'    from peptides, subgroups, groups or taxonomic levels. The format 
#'    of this object is similar to that generated from the 
#'    functions "getsc_specific" and "crumble_taxonomy."
#' 
#' @param target_variable Character indicating the variable that contains 
#'    the elements to be removed. The options are : i) "peptides", ii) "proteins", 
#'    iii) "soubgroups", iv) "groups", v) "sequences", vi) "species", vii) "genus",
#'    viii) "family", ix) "order", x) "class", xi) "phylum" or xii) "superkingdom". 
#'    To select xiii) "samples", it should be indicated the name of ONE column from 
#'    metadata.
#'
#' @param list_elements Atomic vector indicating the elements to be removed. 
#'    For "samples", indicate the element(s) present in the provided variable 
#'    from metadata. For "peptides", "proteins", "subgroups" and "groups" 
#'    provide the X!Tandem nomenclature. For "sequences", provide the 
#'    peptide sequences expressed as aminoacids. For any taxonomic level, provide
#'    the taxonomic entities.
#'
#' @return A list defined as "spectral_count_object" without the elements 
#'    provided in the second argument of the function.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' data(species_fw)
#' 
#' data_selected_samples <- remove_element(fecal_waters, "Methods", c("S_EF", "EF"))
#' 
#' data_selected_peptides <- remove_element(fecal_waters, "peptides", c("pepa3c417", "pepd4664a1"))
#' 
#' data_selected_proteins <- remove_element(species_fw, "proteins", c("a3.a9.a1", "a5.b81.a1"))
#' 
#' data_selected_subgroups <- remove_element(species_fw, "subgroups", c("a3.a9", "b73.a5"))
#' 
#' data_selected_groups <- remove_element(species_fw, "groups", c("a3", "b34", "c231"))
#' 
#' data_selected_sequences <- remove_element(species_fw, "sequences", c("AQLNFGGTIENVVIRDEFPLEK"))
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

remove_element <- 

function(spectral_count_object, target_variable, list_elements){
  
  spectral_count_object <- spectral_count_object
  type_sc_data <- names(spectral_count_object[1])
  
  if (length(spectral_count_object) == 4) {
    
    if (spectral_count_object[[4]] == "spectral_count_object") {
      
      target_variable <- target_variable
      metadata <- spectral_count_object[[2]]
      
      if (target_variable == "peptides" | target_variable == "proteins" | 
          target_variable == "subgroups" | target_variable == "groups" |
          target_variable == "sequences" | target_variable == "species" |
          target_variable == "genus" | target_variable == "family" |
          target_variable == "order" | target_variable == "class" |
          target_variable == "phylum" | target_variable == "superkingdom"){
        
        sc_data <- spectral_count_object[[1]]
        peps_prots <- spectral_count_object[[3]]
        
        #--------------- Removal of list_elements from target_variable
        
        if (target_variable == "peptides") {
          if (all(unique(list_elements) %in% peps_prots$Peptide) == TRUE) {
            new_peps_prots <- peps_prots[-which(peps_prots$Peptide %in% list_elements),]
          }
          else{
            vars <- paste(list_elements, collapse = ", ")
            stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
          }
        }
        else if (target_variable == "proteins") {
          if (all(unique(list_elements) %in% peps_prots$Protein) == TRUE) {
            new_peps_prots <- peps_prots[-which(peps_prots$Protein %in% list_elements),]
          }
          else{
            vars <- paste(list_elements, collapse = ", ")
            stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
          }
        }
        else if (target_variable == "subgroups") {
          if (all(unique(list_elements) %in% peps_prots$Subgroup) == TRUE) {
            new_peps_prots <- peps_prots[-which(peps_prots$Subgroup %in% list_elements),]
          }
          else{
            vars <- paste(list_elements, collapse = ", ")
            stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
          }
        }
        else if (target_variable == "groups") {
          if (all(unique(list_elements) %in% peps_prots$Group) == TRUE) {
            new_peps_prots <- peps_prots[-which(peps_prots$Group %in% list_elements),]
          }
          else{
            vars <- paste(list_elements, collapse = ", ")
            stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
          }
        }
        else if (target_variable == "sequences") {
          if (all(unique(list_elements) %in% peps_prots$Sequence) == TRUE) {
            new_peps_prots <- peps_prots[-which(peps_prots$Sequence %in% list_elements),]
          }
          else{
            vars <- paste(list_elements, collapse = ", ")
            stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
          }
        }
        else if (target_variable == "species") {
          if ("species" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$species) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$species %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'species'")
          }
        }
        else if (target_variable == "genus") {
          if ("genus" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$genus) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$genus %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'genera'")
          }
        }
        else if (target_variable == "family") {
          if ("family" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$family) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$family %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'families'")
          }
        }
        else if (target_variable == "order") {
          if ("order" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$order) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$order %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'order'")
          }
        }
        else if (target_variable == "class") {
          if ("class" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$class) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$class %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'class'")
          }
        }
        else if (target_variable == "phylum") {
          if ("phylum" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$phylum) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$phylum %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'phylum'")
          }
        }
        else if (target_variable == "superkingdom") {
          if ("superkingdom" %in% colnames(peps_prots)){
            if (all(unique(list_elements) %in% peps_prots$superkingdom) == TRUE){
              new_peps_prots <- peps_prots[-which(peps_prots$superkingdom %in% list_elements),]
            }
            else{
              vars <- paste(list_elements, collapse = ", ")
              stop(paste(c("The elements ", vars, " must be present in: ", target_variable), collapse = " "))
            }
          }
          else{
            stop("Verify that the spectral_count_objet is organized by 'superkingdoms'")
          }
        }
        else{
          stop("Verify the provided list of elements")
        }
        
        #--------------- Filtering in spectral counting dataframe with/without taxonomy
        if("species" %in% colnames(peps_prots) | "genus" %in% colnames(peps_prots) | 
           "family" %in% colnames(peps_prots) | "order" %in% colnames(peps_prots) | 
           "class" %in% colnames(peps_prots) | "phylum" %in% colnames(peps_prots) | 
           "superkingdom" %in% colnames(peps_prots)){
          
          if("species" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$species
          }
          else if("genus" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$genus
          }
          else if("family" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$family
          }
          else if("order" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$order
          }
          else if("class" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$class
          }
          else if("phylum" %in% colnames(peps_prots)){
            items_to_keep <- new_peps_prots$phylum
          }
          else{
            items_to_keep <- new_peps_prots$superkingdom
          }
        }
        else if (type_sc_data == "SC_specific_peptides") {
          items_to_keep <- new_peps_prots$Peptide
        }
        else if (type_sc_data == "SC_subgroups") {
          items_to_keep <- new_peps_prots$Subgroup
        }
        else{
          items_to_keep <- new_peps_prots$Group
        }
        
        new_sc_data <- sc_data[rownames(sc_data) %in% items_to_keep,]
        
        #--------------- Output
        # Status message
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
      
      else if (target_variable %in% colnames(metadata) == TRUE) {
        if (all(unique(list_elements) %in% metadata[[target_variable]]) == TRUE) {
          
          #--------------- Logic to remove samples
          # Select SC_samples
          new_metadata <- metadata[-which(metadata[[target_variable]] %in% list_elements),]
          msrun_to_keep <- new_metadata$SC_name
          sc_data <- spectral_count_object[[1]]
          peps_prots <- spectral_count_object[[3]]
          new_sc_data <- sc_data[, colnames(sc_data) %in% msrun_to_keep, drop = FALSE]
          
          # Remove rows with zero spectra
          if(length(unique(msrun_to_keep)) == 1){
            new_sc_data <- new_sc_data[new_sc_data[[msrun_to_keep]] >= 1, , drop = FALSE]
          }
          else{
            elements_to_remove <- NULL
            for (item in 1:dim(new_sc_data)[1]) {
              all_zeros <- all(new_sc_data[item,] == 0)
              if (all_zeros == TRUE) {
                elements_to_remove <- c(elements_to_remove, rownames(new_sc_data[item,]))
              }
            }
            new_sc_data <- new_sc_data[!rownames(new_sc_data) %in% elements_to_remove, , drop = FALSE]
          }
          
          # Remove rows in peps_prots
          if("species" %in% colnames(peps_prots) | "genus" %in% colnames(peps_prots) | 
             "family" %in% colnames(peps_prots) | "order" %in% colnames(peps_prots) | 
             "class" %in% colnames(peps_prots) | "phylum" %in% colnames(peps_prots) | 
             "superkingdom" %in% colnames(peps_prots)){
            if("species" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$species %in% elements_to_remove,]
            }
            else if("genus" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$genus %in% elements_to_remove,]
            }
            else if("family" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$family %in% elements_to_remove,]
            }
            else if("order" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$order %in% elements_to_remove,]
            }
            else if("class" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$class %in% elements_to_remove,]
            }
            else if("phylum" %in% colnames(peps_prots)){
              new_peps_prots <- peps_prots[!peps_prots$phylum %in% elements_to_remove,]
            }
            else{
              new_peps_prots <- peps_prots[!peps_prots$superkingdom %in% elements_to_remove,]
            }
          }
          else if (type_sc_data == "SC_specific_peptides") {
            new_peps_prots <- peps_prots[!peps_prots$Peptide %in% elements_to_remove,]
          }
          else if (type_sc_data == "SC_subgroups") {
            new_peps_prots <- peps_prots[!peps_prots$Subgroup %in% elements_to_remove,]
          }
          else if (type_sc_data == "SC_groups") {
            new_peps_prots <- peps_prots[!peps_prots$Group %in% elements_to_remove,]
          }
          else{
            stop("Invalid object")
          }
          
          # Output
          # Status message
          print("========= Object information ============")
          print(paste("Number of samples: ", length(unique(new_metadata$msrunfile)), sep = ""))
          print(paste("Number of groups: ", length(unique(new_peps_prots$Group)), sep = ""))
          print(paste("Number of subgroups or metaproteins: ", length(unique(new_peps_prots$Subgroup)), sep = ""))
          print(paste("Number of peptides: ", length(unique(new_peps_prots$Peptide)), sep = ""))
          print("=========================================")
          output <- list(new_sc_data, new_metadata, new_peps_prots, "spectral_count_object")
          names(output)[1] <- type_sc_data
          names(output)[2] <- "metadata"
          names(output)[3] <- "peptides_proteins"
          names(output)[4] <- "type_object"
          print("Spectral count object generated")
          return(output)
        }
        else{
          vars <- paste(list_elements, collapse = ", ")
          stop(paste(c("The elements ", vars, " must be present in: '", target_variable, "'. Verify metadata."), collapse = " "))
        }
      }
      else{
        stop("The second argument (target_variable) is invalid, verify this argument")
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