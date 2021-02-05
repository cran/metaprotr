#' add_taxonomy
#'
#' Integrates the database containing the taxonomic classification 
#' of the identified proteins in a "metaproteome_object". 
#' The taxonomic classification is previously obtained by aligment 
#' algorithms and must include seven taxonomic levels assigned 
#' to a given protein: species, genus, family, order, class, phylum 
#' and superkingdom
#' 
#' @param metaproteome_object List defined as "metaproteome_object" 
#'    containing proteins and peptides abundances. The format of this 
#'    object is similar to that generated from the function "load_protspeps".
#' 
#' @param taxonomic_database Dataframe containing the taxonomic information 
#'    for each protein. The first column must contain the same identifiers 
#'    of those present in the column "Accession" from the dataframe 
#'    "peptides_proteins" of the "metaproteome_object". Two additional columns 
#'    have to be present: i) one named "organism" containing the name of the strain
#'    assigned to a given protein; and ii) the other named 
#'    "species.genus.family.order.class.phylum.superkingdom". 
#'    The taxonomic classification can be obtained from a tool of 
#'    sequences aligment and must be ordered as follows: species, genus, 
#'    family, order, class, phylum and superkingdom. The characters inside must 
#'    be concatenated by a comma 
#'    (ex."Streptococcus anginosus,Streptococcus,Streptococcaceae,Lactobacillales,Bacilli,Firmicutes,Bacteria").
#'    An example can be found in this \href{https://zenodo.org/record/3997093}{repository}.
#' 
#' @return A "metaproteome_object", which is a list of six elements with format
#'    similar to that generated from the function "load_protspeps". An additional 
#'    column containing the taxonomic annotation is added to the dataframe 
#'    named "peptides_proteins".
#' 
#' @examples
#' \dontrun{
#' 
#' # Download taxonmical annotation db: https://zenodo.org/record/3997093#.X0UYI6Zb_mE 
#' meta99_full_taxo <- read.csv2("MetaHIT99_best_hit_taxo_complete.tsv", header = TRUE, sep="\t") 
#' 
#' # Files with spectral abundance and proteins list from X!Tandempipeline
#' protein_file <- "your/specific/location/protein_list.txt"
#' peptide_file <- "your/specific/location/peptide_counting.txt"
#' metadata_file <- "your/location/metadata.csv"
#' 
#' metaproteome <- load_protspeps(proteins_file, peptides_file, metadata_file)
#' 
#' metaproteome_taxo <- add_taxonomy(metaproteome, meta99_full_taxo) 
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

add_taxonomy <- 

function(metaproteome_object, taxonomic_database){
  
  metaproteome_object <- metaproteome_object
  
  if (length(metaproteome_object) == 6) {
    if (metaproteome_object[[6]] == "metaproteome_object") {
      peps_prots <- metaproteome_object[[3]]
      
      # Once metaproteome_object validated, check that taxonomic DB was not already integrated
      if ("species.genus.family.order.class.phylum.superkingdom" %in% colnames(peps_prots) == FALSE) {
        print("====================================================")
        taxonomic_database <- taxonomic_database
        print("Assigning taxonomic annotation ...")
        metaproteins_found <- subset(taxonomic_database, taxonomic_database[[1]] %in% peps_prots$Accession)
        percent_found <- (dim(metaproteins_found)[1] / dim(peps_prots[which(!duplicated(peps_prots$Protein)),])[1]) * 100
        print(paste(dim(metaproteins_found)[1], "proteins with taxonomy", sep = " "))
        print(paste(round(percent_found, 2), "% of assigned proteins from the dataset", sep = " "))
        new_peps_prots <- merge(
          x = peps_prots,
          y = metaproteins_found[, c(
            colnames(metaproteins_found)[1],
            "organism",
            "species.genus.family.order.class.phylum.superkingdom"
          )],
          by.x = "Accession",
          by.y = colnames(metaproteins_found)[1],
          all = TRUE
        )
        new_peps_prots$organism <- as.character(new_peps_prots$organism)
        new_peps_prots$species.genus.family.order.class.phylum.superkingdom <-
          as.character(new_peps_prots$species.genus.family.order.class.phylum.superkingdom)
        new_peps_prots[["species.genus.family.order.class.phylum.superkingdom"]][is.na(new_peps_prots[["species.genus.family.order.class.phylum.superkingdom"]])] <- "unassigned,unassigned,unassigned,unassigned,unassigned,unassigned,unassigned"
        new_peps_prots[["organism"]][is.na(new_peps_prots[["organism"]])] <- "unassigned"
        
        proteins_info <- metaproteome_object[[1]]
        pepdata <- metaproteome_object[[2]]
        metadata <- metaproteome_object[[4]]
        SCpep <- metaproteome_object[[5]]
        type_object <- metaproteome_object[[6]]

        # Status message
        print("Creating metaproteome object")
        print("========= Object information ============")
        print(paste("Number of samples: ", length(unique(metadata$msrunfile)), sep=""))
        print(paste("Number of groups: ", length(unique(new_peps_prots$Group)), sep=""))
        print(paste("Number of subgroups or metaproteins: ", length(unique(new_peps_prots$Subgroup)), sep=""))
        print(paste("Number of peptides: ", length(unique(new_peps_prots$Peptide)), sep=""))
        print("=========================================")
        
        #  OUTPUT DATA  
        output <-
          list(proteins_info,
               pepdata,
               new_peps_prots,
               metadata,
               SCpep,
               type_object)
        names(output)[1] <- "proteins_data"
        names(output)[2] <- "peptides_data"
        names(output)[3] <- "peptides_proteins"
        names(output)[4] <- "metadata"
        names(output)[5] <- "spectral_counting_peptides"
        names(output)[6] <- "type_object"
        print("metaproteome object generated")
        return(output)
        print("====================================================")
      }
      else{
        stop("The taxonomic database was already loaded to the metaproteome object")
      }
    }
    else {
      stop("Invalid metaproteome_object object")
    }
  }
  else{
    stop("Invalid metaproteome_object object")
  }
}