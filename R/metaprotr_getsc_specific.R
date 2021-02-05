#' getsc_specific
#'
#' Returns the abundances, expressed as spectral counts (SC), of the different 
#' peptides, subgroups (also referred as metaprotein) or groups within the 
#' samples of the experiment. The abundance corresponds to the sum of SC of the 
#' specific peptides present in a given subgroup or group. See 
#' \href{http://pappso.inrae.fr/bioinfo/xtandempipeline/}{X!TandemPipeline} 
#' for more details concerning the grouping algorithm.
#' 
#' @param metaproteome_object List defined as "metaproteome_object" with 
#'    dataframes having a similar format to that generated from 
#'    "load_protspeps" function. It contains metaproteomics 
#'    data such as peptide and protein abundances and sample 
#'    information.
#' 
#' @param type_SCspecific Character indicating the type of data 
#'    to be returned. The possible options are: i) "sc_specific_peptides" 
#'    for peptides, or ii) "sc_groups" for groups, or iii) "sc_subgroups" 
#'    for subgroups. For more details of the identification algorithm 
#'    (peptides, subgroups, groups) check 
#'    \href{http://pappso.inrae.fr/bioinfo/xtandempipeline/}{X!TandemPipeline}.
#' 
#' @return A list of four elements defined as "spectral_count_object". The first 
#'    element is a dataframe organized in function of peptides OR subgroups 
#'    (also referred as metaproteins) OR groups. The entities of this dataframe 
#'    have their abundance expressed as SC from specific peptides. The second 
#'    element is a dataframe of metadata containing the experiment information. 
#'    The third element is a dataframe containing the information of peptides 
#'    with their associated proteins. The fourth element is a character 
#'    indicating the type of object generated.
#' 
#' @examples
#' \dontrun{
#' 
#' # From a given "metaproteome_object" add the taxonomic classification
#' 
#' metaproteome <- load_protspeps(proteins_file, peptides_file, metadata_file)
#' metaproteome_taxo <- add_taxonomy(metaproteome, meta99_full_taxo) 
#' 
#' # Organize proteomics data by peptides OR subgroups OR groups
#' 
#' SC_specific_peptides <- getsc_specific(metaproteome_taxo, 'sc_specific_peptides') 
#' SC_specific_groups <- getsc_specific(metaproteome_taxo, 'sc_groups') 
#' SC_specific_subgroups <- getsc_specific(metaproteome_taxo, 'sc_subgroups') 
#' 
#' }
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

getsc_specific <- 

function(metaproteome_object, type_SCspecific){
  metaproteome_object <- metaproteome_object
  
  if ((length(metaproteome_object) == 6) && (metaproteome_object[[6]] == "metaproteome_object")) {
    peps_prots <- metaproteome_object[[3]]
    metadata <- metaproteome_object[[4]]
    SCpeps <- metaproteome_object[[5]]
    spec <- peps_prots$Peptide[ave(peps_prots$Peptide, peps_prots$Peptide, FUN = length) == 1] #specific peptides
    SCpeps_specific <- SCpeps[spec, ]
    peps_prots_specific <- peps_prots[peps_prots$Peptide %in% spec,]
    
    if (type_SCspecific == "sc_specific_peptides") {
      nb_elements <- dim(SCpeps_specific)[1]
      output <-
        list(SCpeps_specific,
             metadata,
             peps_prots_specific,
             "spectral_count_object")
      names(output)[1] <- "SC_specific_peptides"
      names(output)[2] <- "metadata"
      names(output)[3] <- "peptides_proteins"
      names(output)[4] <- "type_object"
      print(paste("Spectral data with ", as.character(nb_elements), " specific peptides", sep = ""))
      return(output)
    }
    
    else if (type_SCspecific == "sc_subgroups") {
      SCpeps_specific_subgroups <-
        merge(SCpeps_specific,
              peps_prots[, c("Subgroup", "Peptide")],
              by.x = "row.names",
              by.y = "Peptide")
      SCpeps_specific_subgroups <- SCpeps_specific_subgroups[, 2:ncol(SCpeps_specific_subgroups)]
      SCpeps_specific_subgroups_combine <- aggregate(. ~ Subgroup, SCpeps_specific_subgroups, sum)
      rownames(SCpeps_specific_subgroups_combine) <- SCpeps_specific_subgroups_combine[, "Subgroup"]
      SCsubgroups <- SCpeps_specific_subgroups_combine[2:ncol(SCpeps_specific_subgroups_combine)]
      nb_elements <- dim(SCsubgroups)[1]
      output <-
        list(SCsubgroups,
             metadata,
             peps_prots_specific,
             "spectral_count_object")
      names(output)[1] <- "SC_subgroups"
      names(output)[2] <- "metadata"
      names(output)[3] <- "peptides_proteins"
      names(output)[4] <- "type_object"
      print(paste("Spectral data with ", as.character(nb_elements), " subgroups or metaproteins", sep = ""))
      return(output)
    }
    
    else if (type_SCspecific == "sc_groups") {
      groups_peps <- peps_prots[, c("Group", "Peptide")]
      groups_peps <- groups_peps[!duplicated(groups_peps), ]
      SCpeps_groups <- merge(SCpeps, groups_peps, by.x = "row.names", by.y = "Peptide")
      SCpeps_groups <- SCpeps_groups[, 2:ncol(SCpeps_groups)]
      SCpeps_groups_combined <- aggregate(. ~ Group, SCpeps_groups, sum)
      rownames(SCpeps_groups_combined) <- SCpeps_groups_combined[, "Group"]
      SCgroups <- SCpeps_groups_combined[2:ncol(SCpeps_groups_combined)]
      nb_elements <- dim(SCgroups)[1]
      output <-
        list(SCgroups,
             metadata,
             peps_prots_specific,
             "spectral_count_object")
      names(output)[1] <- "SC_groups"
      names(output)[2] <- "metadata"
      names(output)[3] <- "peptides_proteins"
      names(output)[4] <- "type_object"
      print(paste("Spectral data with ", as.character(nb_elements), " groups", sep = ""))
      return(output)
    }
    
    else{
      stop(
        "Select ONE of the following options to export: 'sc_specific_peptides', 'sc_subgroups' or 'sc_groups'"
      )
    }
  }
  
  else{
    stop("Invalid metaproteome object")
  } 
}

