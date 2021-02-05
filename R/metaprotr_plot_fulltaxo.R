#' plot_fulltaxo
#' 
#' Provides the number of taxonomic entities per sample in the different 
#' taxonomic levels. The taxonomic levels are: species, genus, family, order, 
#' class, phylum and superkingdom.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with format similar to that generated with 
#'    the function "getsc_specific".  This object contains abundances 
#'    expressed as spectral counts from peptides, subgroups (metaproteins) 
#'    or groups. Taxonomy must have previously been added with the function 
#'    "add_taxonomy".
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a pdf and a csv file in the workstation of 
#'    the user.
#' 
#' @return Bar plots (pdf) and csv file with the number of taxonomic species, 
#'    genus, family, order, class, phylum and superkingdom per sample. 
#'    An additional csv file is generated providing the rate of assignment. 
#'    The rate of assigment corresponds to the ratio between the number of 
#'    the most frequent annotation ("species", "genus", "family", "order", 
#'    "class", "phylum" or "superkingdom") and the total number of elements within 
#'    each level of the spectral category under study (subgroup or group).
#'    The csv file is generated only when the "spectral_count_object" is 
#'    organized by subgroup or by group.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' plot_fulltaxo(fecal_waters)
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

plot_fulltaxo <-

function(spectral_count_object, force = FALSE){
  
  counts <- genus <- Group <- Peptide <- phylum <- SampleID <- species <- Subgroup <- superkingdom <- taxonomy <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <- select.list(c("yes", "no"), 
                            title = "Do you allow to create a pdf and a csv file with taxonomic information?")
    if (response == "yes") {
      if (length(spectral_count_object) == 4) {
        if(spectral_count_object[[4]] == "spectral_count_object") {
          
          sc_data <- spectral_count_object[[1]]
          type_sc_data <- names(spectral_count_object)[1]
          metadata <- spectral_count_object[[2]]
          pepsprots <- spectral_count_object[[3]]
          
          if("species.genus.family.order.class.phylum.superkingdom" %in% colnames(pepsprots)){
            print("============== Creating object ==============")
            # Split taxonomy
            taxo <- colsplit(pepsprots$species.genus.family.order.class.phylum.superkingdom, 
                             ",", c("species", "genus","family", "order", "class", "phylum", "superkingdom"))
            new_pepsprots <- cbind(pepsprots, taxo)
            new_pepsprots$species <- as.factor(new_pepsprots$species)
            new_pepsprots$genus <- as.factor(new_pepsprots$genus)
            new_pepsprots$family <- as.factor(new_pepsprots$family)
            new_pepsprots$order <- as.factor(new_pepsprots$order)
            new_pepsprots$class <- as.factor(new_pepsprots$class)
            new_pepsprots$phylum <- as.factor(new_pepsprots$phylum)
            new_pepsprots$superkingdom <- as.factor(new_pepsprots$superkingdom)
            
            # Count of hits of the most frequent factor per taxonomic level
            # in order to know the % of assignment per level
            wrapp_counts <- NULL
            print("Calculating counts ...")
            print("This step will take some time be patient ...")
            
            if(type_sc_data == "SC_groups"){
              
              # Logic for groups
              spectral_origin <- "Group"
              elements <- unique(levels(as.factor(new_pepsprots$Group)))
              for(i in 1:length(unique(new_pepsprots$Group))){
                grouped_data <- new_pepsprots[new_pepsprots$Group == elements[i], ]
                grouped_data <- droplevels(grouped_data)
                species_max_hits <- max(summary(grouped_data$species))
                grouped_data <- droplevels(grouped_data)
                genus_max_hits <- max(summary(grouped_data$genus))
                grouped_data <- droplevels(grouped_data)
                family_max_hits <- max(summary(grouped_data$family))
                grouped_data <- droplevels(grouped_data)
                order_max_hits <- max(summary(grouped_data$order))
                grouped_data <- droplevels(grouped_data)
                class_max_hits <- max(summary(grouped_data$class))
                grouped_data <- droplevels(grouped_data)
                phylum_max_hits <- max(summary(grouped_data$phylum))
                grouped_data <- droplevels(grouped_data)
                superkingdom_max_hits <- max(summary(grouped_data$superkingdom))
                grouped_data <- droplevels(grouped_data)
                wrapp_counts <- rbind.data.frame(wrapp_counts, cbind.data.frame(
                  Group = elements[i], 
                  hits_species = species_max_hits, 
                  hits_genus = genus_max_hits,
                  hits_family = family_max_hits,
                  hits_order = order_max_hits,
                  hits_class = class_max_hits,
                  hits_phylum = phylum_max_hits,
                  hits_superkingdom = superkingdom_max_hits))
              }
              # Taxonomic levels
              taxo_assign_max <- new_pepsprots %>% group_by(Group)  %>%  summarize(
                species = names(which.max(table(species))),
                genus = names(which.max(table(genus))),
                family = names(which.max(table(family))),
                order = names(which.max(table(order))),
                class = names(which.max(table(class))),
                phylum = names(which.max(table(phylum))),
                superkingdom = names(which.max(table(superkingdom))),
                total_nb = table(Group))
              wrapp_counts$Group <- as.character(wrapp_counts$Group)
              # Export of rate of assignment by groups
              wrapp_counts <- merge(wrapp_counts, taxo_assign_max, by=spectral_origin)
              wrapp_counts[, c(2:8)] <- round(wrapp_counts[, c(2:8)] / wrapp_counts$total_nb, 3)
              filename <- paste(spectral_origin, "_taxonomy_assignment.csv", sep = "")
              write.table(wrapp_counts, filename, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
              print(paste("CSV file ", filename, "generated", sep = " "))
            }
            else if(type_sc_data == "SC_subgroups"){
              # Logic for subgroups
              spectral_origin <- "Subgroup"
              elements <- unique(levels(as.factor(new_pepsprots$Subgroup)))
              for(i in 1:length(unique(new_pepsprots$Subgroup))){
                grouped_data <- new_pepsprots[new_pepsprots$Subgroup == elements[i], ]
                grouped_data <- droplevels(grouped_data)
                species_max_hits <- max(summary(grouped_data$species))
                grouped_data <- droplevels(grouped_data)
                genus_max_hits <- max(summary(grouped_data$genus))
                grouped_data <- droplevels(grouped_data)
                family_max_hits <- max(summary(grouped_data$family))
                grouped_data <- droplevels(grouped_data)
                order_max_hits <- max(summary(grouped_data$order))
                grouped_data <- droplevels(grouped_data)
                class_max_hits <- max(summary(grouped_data$class))
                grouped_data <- droplevels(grouped_data)
                phylum_max_hits <- max(summary(grouped_data$phylum))
                grouped_data <- droplevels(grouped_data)
                superkingdom_max_hits <- max(summary(grouped_data$superkingdom))
                grouped_data <- droplevels(grouped_data)
                wrapp_counts <- rbind.data.frame(wrapp_counts, cbind.data.frame(
                  Subgroup = elements[i], 
                  hits_species = species_max_hits, 
                  hits_genus = genus_max_hits,
                  hits_family = family_max_hits,
                  hits_order = order_max_hits,
                  hits_class = class_max_hits,
                  hits_phylum = phylum_max_hits,
                  hits_superkingdom = superkingdom_max_hits))
              }
              # Taxonomic levels
              taxo_assign_max <- new_pepsprots %>% group_by(Subgroup)  %>%  summarize(
                species = names(which.max(table(species))),
                genus = names(which.max(table(genus))),
                family = names(which.max(table(family))),
                order = names(which.max(table(order))),
                class = names(which.max(table(class))),
                phylum = names(which.max(table(phylum))),
                superkingdom = names(which.max(table(superkingdom))),
                total_nb = table(Subgroup))
              wrapp_counts$Subgroup <- as.character(wrapp_counts$Subgroup)
              # Export of rate of assignment by subgroups
              wrapp_counts <- merge(wrapp_counts, taxo_assign_max, by=spectral_origin)
              wrapp_counts[, c(2:8)] <- round(wrapp_counts[, c(2:8)] / wrapp_counts$total_nb, 3)
              filename <- paste(spectral_origin, "_taxonomy_assignment.csv", sep = "")
              write.table(wrapp_counts, filename, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
              print(paste("CSV file ", filename, "generated", sep = " "))
            }
            else{
              
              # Logic for peptides
              spectral_origin <- "Peptide"
              # Taxonomic levels
              taxo_assign_max <- new_pepsprots %>% group_by(Peptide)  %>%  summarize(
                species = names(which.max(table(species))),
                genus = names(which.max(table(genus))),
                family = names(which.max(table(family))),
                order = names(which.max(table(order))),
                class = names(which.max(table(class))),
                phylum = names(which.max(table(phylum))),
                superkingdom = names(which.max(table(superkingdom))),
                total_nb = table(Peptide))
            }
            
            # Create dataframe with counts per taxonomic level
            species_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "species")])
            species_tab <- apply(species_tab[, 1:length(colnames(sc_data))], 2, tapply,
                                 species_tab$species, sum)
            species_counts <- colSums(species_tab > 0)
            species_counts <- as.data.frame(enframe(species_counts, name = "SC_name", value = "counts"))
            species_counts$taxonomy <- "1_species"
            genus_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "genus")])
            genus_tab <- apply(genus_tab[, 1:length(colnames(sc_data))], 2, tapply,
                               genus_tab$genus, sum)
            genus_counts <- colSums(genus_tab > 0)
            genus_counts <- as.data.frame(enframe(genus_counts, name = "SC_name", value = "counts"))
            genus_counts$taxonomy <- "2_genus"
            family_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "family")])
            family_tab <- apply(family_tab[, 1:length(colnames(sc_data))], 2, tapply,
                                family_tab$family, sum)
            family_counts <- colSums(family_tab > 0)
            family_counts <- as.data.frame(enframe(family_counts, name = "SC_name", value = "counts"))
            family_counts$taxonomy <- "3_family"
            order_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "order")])
            order_tab <- apply(order_tab[, 1:length(colnames(sc_data))], 2, tapply,
                               order_tab$order, sum)
            order_counts <- colSums(order_tab > 0)
            order_counts <- as.data.frame(enframe(order_counts, name = "SC_name", value = "counts"))
            order_counts$taxonomy <- "4_order"
            class_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "class")])
            class_tab <- apply(class_tab[, 1:length(colnames(sc_data))], 2, tapply,
                               class_tab$class, sum)
            class_counts <- colSums(class_tab > 0)
            class_counts <- as.data.frame(enframe(class_counts, name = "SC_name", value = "counts"))
            class_counts$taxonomy <- "5_class"
            phylum_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "phylum")])
            phylum_tab <- apply(phylum_tab[, 1:length(colnames(sc_data))], 2, tapply,
                                phylum_tab$phylum, sum)
            phylum_counts <- colSums(phylum_tab > 0)
            phylum_counts <- as.data.frame(enframe(phylum_counts, name = "SC_name", value = "counts"))
            phylum_counts$taxonomy <- "6_phylum"
            superkingdom_tab <- cbind(sc_data, taxo_assign_max[,c(spectral_origin, "superkingdom")])
            superkingdom_tab <- apply(superkingdom_tab[, 1:length(colnames(sc_data))], 2, tapply,
                                      superkingdom_tab$superkingdom, sum)
            superkingdom_counts <- colSums(superkingdom_tab > 0)
            superkingdom_counts <- as.data.frame(enframe(superkingdom_counts, name = "SC_name", value = "counts"))
            superkingdom_counts$taxonomy <- "7_superkingdom"
            
            taxo_counts <- rbind(species_counts, genus_counts, family_counts, 
                                 order_counts, class_counts, phylum_counts, superkingdom_counts)
            taxo_counts <- merge(taxo_counts, metadata[, c("SC_name", "SampleID")], by="SC_name")
            taxo_counts$SC_name <- as.factor(taxo_counts$SC_name)
            taxo_counts$SampleID <- as.factor(taxo_counts$SampleID)
            taxo_counts$taxonomy <- as.factor(taxo_counts$taxonomy)
            
            # Export counts
            filename <- paste(spectral_origin, "_taxonomy_counts.csv", sep = "")
            write.table(taxo_counts, filename, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
            print(paste("CSV file ", filename, "generated", sep = " "))
            
            # To keep graphics settings
            old_par <- par(no.readonly = TRUE)
            on.exit(suppressWarnings(par(old_par)))
            
            ylimits <- max(taxo_counts$counts) + 30
            taxo_plot <- ggplot(data = taxo_counts, aes(x = taxonomy, y = counts, fill = SampleID)) + 
              geom_bar(stat="identity", position = position_dodge()) +
              geom_text(aes(label = counts, angle = 90), position = position_dodge(width = 0.9), vjust = 0.5, hjust = -0.2, size = 3.5) +
              labs(title = paste("Number of taxons per taxonomic level and per sample ( ", spectral_origin, " matrix)", sep = ""), 
                   x = "", 
                   y = "Elements") +
              ylim(-1, ylimits) +
              scale_fill_viridis_d() +
              theme_classic(base_size = 19)
            filename <- paste(spectral_origin, "_taxonomy_count.pdf", sep = "")
            pdf(filename, width = 16, height = 10)
            print(taxo_plot)
            dev.off()
            print(paste("The figure", filename, "was generated", sep = " "))
          }
          else{
            stop("Invalid spectral_count_object. Taxonomy must be added with the function add_taxonomy")
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
    else{
      stop("No file was created")
    }
  }
}