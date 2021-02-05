#' crumble_taxonomy
#'
#' Generates a list of four elements defined as "spectral_count_object"
#' containing taxonomic classification. The first element is a dataset 
#' that contains the spectral counts abundance organized by a provided 
#' taxonomic level. The possible taxonomic levels are: species, genus, 
#' family, order, class, phylum or superkingdom.
#' 
#' @param spectral_count_object List defined as "spectral_count_object"
#'    containing dataframes with abundance expressed as spectral counts
#'    and organized by peptides, subgroups or groups. The format of 
#'    this object is similar to that generated with the function 
#'    "getsc_specific". Taxonomy must be added previously 
#'    with "add_taxonomy" function.
#'
#' @param taxonomic_level Character indicating the taxonomic level to which
#'    the spectral abundance will be arranged in the samples of the 
#'    "spectral_count_object". The possible options are: "species", "genus", 
#'    "family", "order", "class", "phylum" or "superkingdom".
#'
#' @param filter_rate Numeric value between 0 and 1 that indicates the 
#'    minimal rate of consensual annotation desired by the user within each 
#'    level of the spectral category (subgroup or group). This rate is defined 
#'    as the ratio between the number of the most frequent annotation entity 
#'    ("species", "genus", "family", "order", "class", "phylum" or "superkingdom") 
#'    divided by the total number of entities within each level of the spectral 
#'    category under study (subgroup or group). The default value is set to 1, 
#'    if 100 \% of consensus is desired.
#'
#' @return A list of four elements defined as "spectral_count_object", the first 
#'    element is a dataframe with abundance expressed as spectral counts of entities 
#'    (peptides, subgroups or groups) organized by the provided taxonomic level. 
#'    The second element is a dataframe that contains the experiment information. 
#'    The third element is a dataframe containing the information of peptides with 
#'    their associated proteins. And the fourth element is a character indicating 
#'    the type of object generated. 
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' 
#' superkingdom_fecalwaters <- crumble_taxonomy(fecal_waters, "superkingdom")
#' 
#' phylum_fecalwaters <- crumble_taxonomy(fecal_waters, "phylum")
#' 
#' class_fecalwaters <- crumble_taxonomy(fecal_waters, "class")
#' 
#' order_fecalwaters <- crumble_taxonomy(fecal_waters, "order")
#' 
#' family_fecalwaters <- crumble_taxonomy(fecal_waters, "family")
#' 
#' genus_fecalwaters <- crumble_taxonomy(fecal_waters, "genus")
#' 
#' species_fecalwaters <- crumble_taxonomy(fecal_waters, "species")
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

crumble_taxonomy <-

function(spectral_count_object, taxonomic_level, filter_rate = 1){
  
  genus <- Group <- Peptide <- phylum <- species <- Subgroup <- superkingdom <- NULL
  spectral_count_object <- spectral_count_object
  
  if(length(spectral_count_object) == 4){
    if(spectral_count_object[[4]] == "spectral_count_object"){
      if(taxonomic_level == "species" | taxonomic_level == "genus" |
         taxonomic_level == "family" | taxonomic_level == "order" | 
         taxonomic_level == "class" | taxonomic_level == "phylum" | 
         taxonomic_level == "superkingdom"){
      
        sc_data <- spectral_count_object[[1]]
        type_sc_data <- names(spectral_count_object)[1]
        metadata <- spectral_count_object[[2]]
        pepsprots <- spectral_count_object[[3]]
        
        if("species.genus.family.order.class.phylum.superkingdom" %in% colnames(pepsprots)){
          if("species" %in% colnames(pepsprots) | "genus" %in% colnames(pepsprots) | 
             "family" %in% colnames(pepsprots) | "order" %in% colnames(pepsprots) | 
             "class" %in% colnames(pepsprots) | "phylum" %in% colnames(pepsprots) | 
             "superkingdom" %in% colnames(pepsprots)){
            stop("Invalid spectral_count_object. The object was already subjected to the function crumble_taxonomy")
          }
          else{
            if(filter_rate < 0 | filter_rate > 1){
              stop("The third argument must a number between 0 and 1")
            }
            else{
              print("Creating object ...")
              #-------------- Split taxonomy
              taxo <- colsplit(pepsprots$species.genus.family.order.class.phylum.superkingdom,
                               ",", c("species", "genus","family", "order", "class", "phylum", "superkingdom"))
              new_pepsprots <- cbind(pepsprots, taxo[, c(taxonomic_level)])
              colnames(new_pepsprots)[length(colnames(new_pepsprots))] <- taxonomic_level
              
              if(type_sc_data == "SC_groups"){
                #-------------- logic for groups
                spectral_origin <- "Group"
                if(taxonomic_level == "species"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    species = names(which.max(table(species))), 
                    total_nb = table(Group))
                }
                else if(taxonomic_level == "genus"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    genus = names(which.max(table(genus))), 
                    total_nb = table(Group))
                }
                else if(taxonomic_level == "family"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    family = names(which.max(table(family))), 
                    total_nb = table(Group))
                }
                else if(taxonomic_level == "order"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    order = names(which.max(table(order))), 
                    total_nb = table(Group))
                }
                else if(taxonomic_level == "class"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    class = names(which.max(table(class))), 
                    total_nb = table(Group))
                }
                else if(taxonomic_level == "phylum"){
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    phylum = names(which.max(table(phylum))), 
                    total_nb = table(Group))
                }
                else{
                  taxo_assign <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                    superkingdom = names(which.max(table(superkingdom))), 
                    total_nb = table(Group))
                }
              }
              else if(type_sc_data == "SC_subgroups"){
                #-------------- logic for subgroups
                spectral_origin <- "Subgroup"
                if(taxonomic_level == "species"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    species = names(which.max(table(species))), 
                    total_nb = table(Subgroup))
                }
                else if(taxonomic_level == "genus"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    genus = names(which.max(table(genus))), 
                    total_nb = table(Subgroup))
                }
                else if(taxonomic_level == "family"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    family = names(which.max(table(family))), 
                    total_nb = table(Subgroup))
                }
                else if(taxonomic_level == "order"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    order = names(which.max(table(order))), 
                    total_nb = table(Subgroup))
                }
                else if(taxonomic_level == "class"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    class = names(which.max(table(class))), 
                    total_nb = table(Subgroup))
                }
                else if(taxonomic_level == "phylum"){
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    phylum = names(which.max(table(phylum))), 
                    total_nb = table(Subgroup))
                }
                else{
                  taxo_assign <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                    superkingdom = names(which.max(table(superkingdom))), 
                    total_nb = table(Subgroup))
                }
              }
              else{
                spectral_origin <- "Peptide"
                if(taxonomic_level == "species"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    species = names(which.max(table(species))), 
                    total_nb = table(Peptide))
                }
                else if(taxonomic_level == "genus"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    genus = names(which.max(table(genus))), 
                    total_nb = table(Peptide))
                }
                else if(taxonomic_level == "family"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    family = names(which.max(table(family))), 
                    total_nb = table(Peptide))
                }
                else if(taxonomic_level == "order"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    order = names(which.max(table(order))), 
                    total_nb = table(Peptide))
                }
                else if(taxonomic_level == "class"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    class = names(which.max(table(class))), 
                    total_nb = table(Peptide))
                }
                else if(taxonomic_level == "phylum"){
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    phylum = names(which.max(table(phylum))), 
                    total_nb = table(Peptide))
                }
                else{
                  taxo_assign <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                    superkingdom = names(which.max(table(superkingdom))), 
                    total_nb = table(Peptide))
                }
              }
              taxo_assign <- as.data.frame(taxo_assign)
              
              #-------------- To apply filter of % assignment to subgroup or groups of metaproteins
              if(spectral_origin == "Group" | spectral_origin == "Subgroup"){
                wrapp_counts <- NULL
                elements <- unique(levels(as.factor(new_pepsprots[[spectral_origin]])))
                print("This step will take some time be patient ...")
                for(i in 1:length(unique(new_pepsprots[[spectral_origin]]))){
                  grouped_data <- new_pepsprots[new_pepsprots[[spectral_origin]] == elements[i], ]
                  grouped_data <- droplevels(grouped_data)
                  max_hits <- max(summary(as.factor(grouped_data[[taxonomic_level]])))
                  wrapp_counts <- rbind.data.frame(wrapp_counts, cbind.data.frame(
                    origin = elements[i], 
                    max_hits = max_hits))
                }
                names(wrapp_counts)[1] <- spectral_origin
                taxo_assign <- merge(taxo_assign, wrapp_counts, by = spectral_origin)
                taxo_assign$assignment <- as.numeric(taxo_assign$max_hits / taxo_assign$total_nb)
              }
              else{
                # Case for peptides
                taxo_assign$assignment <- 1
              }
              
              #-------------- Create dataframe with counts per taxonomic level
              taxo_tab <- merge(
                sc_data, 
                taxo_assign[,c(spectral_origin, taxonomic_level, "assignment")], 
                by.x = "row.names", 
                by.y = spectral_origin)
              taxo_tab <- taxo_tab[taxo_tab$assignment >= filter_rate, ]
              taxo_tab[[taxonomic_level]] <- as.character(taxo_tab[[taxonomic_level]])
              taxo_tab <- apply(taxo_tab[, 2:(length(colnames(taxo_tab))-2)],
                                2, tapply, taxo_tab[[taxonomic_level]], sum)
              taxo_tab <- as.data.frame(taxo_tab)
              new_pepsprots[[taxonomic_level]] <- as.character(new_pepsprots[[taxonomic_level]])
              
              #-------------- Status message
              print("========= Object information ============")
              print(paste("Number of samples: ", length(unique(metadata$msrunfile)), sep = ""))
              print(paste("Taxonomic level: ", taxonomic_level, sep = ""))
              print(paste("Number of elements in level: ", dim(taxo_tab)[1], sep = ""))
              print(paste("Spectral counting origin: ", spectral_origin, sep = ""))
              print("=========================================")
              
              #-------------- OUTPUT DATA 
              output <- list(taxo_tab, metadata, new_pepsprots, "spectral_count_object")
              names(output)[1] <- type_sc_data
              names(output)[2] <- "metadata"
              names(output)[3] <- "peptides_proteins"
              names(output)[4] <- "type_object"
              print(paste("spectral_count_object based on <", taxonomic_level ,"> was created", sep = " "))
              return(output)
            }
          }
        }
        else{
          stop("Invalid spectral_count_object. Taxonomy must be added with the function add_taxonomy")
        }
      }
      else{
        stop("Invalid taxonomic level. The second argument must be : 'species', 'genus', 'family', 'order', 'class', 'phylum' or 'superkingdom'")
      }
    }
    else{
      stop("Invalid spectral_count_object")
    }
  }
  else{
    stop("Invalid spectral_count_object")
  }
}