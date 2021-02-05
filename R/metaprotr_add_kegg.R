#' add_kegg
#' 
#' Integrates a database containing the functional annotation of the identified 
#' metaproteins into a list defined as "spectral_count_object". The proteins 
#' from the “spectral_count_object” must contain taxonomic information. 
#' The functional annotation was obtained from the Kyoto Encyclopedia of 
#' Genes and Genomes (KEGG) Orthology database. This database contains 
#' the molecular functions represented in terms of functional orthologs (KO terms).
#' Check \href{https://www.genome.jp/kegg/}{KEGG} for more details.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing the abundance of the elements (groups, subgroups 
#'    or peptides) expressed as spectral counts and organized by 
#'    taxonomic levels. The format of this object is similar to 
#'    that generated from the function "crumble_taxonomy".
#' 
#' @param annotation_db Dataframe containing the functional annotation 
#'    of the proteins. This dataframe must contain two variables: 
#'    i) "gene_name": indicating the same protein names to those present in 
#'    the variable  "Accession" from the "peptides_proteins", third dataframe in the 
#'    list defined as "spectral_count_object"; and, ii) "ko": indicating the 
#'    KEGG Orthology code assigned to a given protein. An example can be found in
#'    this \href{https://zenodo.org/record/3997093}{repository}.
#' 
#' @param taxonomic_db Dataframe containing the taxonomic information 
#'    for each protein. The first column must contain the same identifiers 
#'    of those present in the column "Accession" from the dataframe 
#'    "peptides_proteins" of the "metaproteome_object". Two additional columns 
#'    have to be present: i) one named "organism" containing the name of the strain
#'    assigned to a given protein; and ii) the other named 
#'    "species.genus.family.order.class.phylum.superkingdom". 
#'    The taxonomic classification can be obtained from a tool of 
#'    sequences aligment and must be ordered as follows: species, genus, 
#'    family, order, class, phylum and superkingdom. The characters inside must 
#'    be concatenated by a comma without spaces 
#'    (ex. "Streptococcus anginosus,Streptococcus,Streptococcaceae,Lactobacillales,Bacilli,Firmicutes,Bacteria").
#'    An example can be found in this \href{https://zenodo.org/record/3997093}{repository}.
#' 
#' @param metaproteome_origin List defined as "metaproteome_object" generated
#'    from the function 'load_protspeps'.
#' 
#' @param protein_file Character indicating the location of a txt 
#'    file containing the list of proteins generated in 
#'    \href{http://pappso.inrae.fr/bioinfo/xtandempipeline/}{X!TandemPipeline} 
#'    using an adapted iterative approach described by 
#'    \href{https://www.theses.fr/2019SORUS043}{Bassignani, 2019}. 
#'    Separation between columns should be indicated by tabulation. 
#'    For more details regarding data input check 
#'    \href{https://forgemia.inra.fr/pappso/metaprotr#data-inputs}{format examples}.
#' 
#' @param peptide_file Character indicating the location of a txt 
#'    file containing peptides abundances expressed as spectral counts.
#'    This file is generated from 
#'    \href{http://pappso.inrae.fr/bioinfo/xtandempipeline/}{X!TandemPipeline} 
#'    using an adapted iterative approach described by 
#'    \href{https://www.theses.fr/2019SORUS043}{Bassignani, 2019}. 
#'    Separation between columns should be indicated by tabulation. 
#'    For more details regarding data input check 
#'    \href{https://forgemia.inra.fr/pappso/metaprotr#data-inputs}{format examples}.
#' 
#' @param text_to_filter Character containig a part of text to be searched in 
#'    the "Description" of the protein file. All the elements containing this 
#'    character will be removed. The default value was set to "HUMAN".
#' 
#' @param taxonomic_levels_allowed Numeric value indicating the maximal
#'    number of taxonomic levels allowed per spectral group or subgroup (in function
#'    of the type of spectral data). The default value is set to 1.
#' 
#' @return A list defined as "spectral_count_object" with the functional 
#'    annotation added to the identified proteins. A new column is added to 
#'    the dataframe "peptides_proteins". Two quality control plot are also generated,
#'    one with the number of taxonomic entities per spectral level and another with
#'    the number of KO terms per spectral level. 
#' 
#' @examples
#' \dontrun{
#' 
#' # Download functional and taxonmical annotation db: https://zenodo.org/record/3997093#.X0UYI6Zb_mE 
#' meta99_full_taxo <- read.csv2("full_taxonomy_MetaHIT99.tsv", header= TRUE, sep="\t")
#' kegg_db <- read.csv2("hs_9_9_igc_vs_kegg89.table", header = TRUE, sep = "\t")
#' 
#' # Files with spectral abundance and proteins list from X!Tandempipeline
#' protein_file <- "your/specific/location/protein_list.txt"
#' peptide_file <- "your/specific/location/peptide_counting.txt"
#' metadata_file <- "your/location/metadata.csv"
#' 
#' metaproteome_origin <- load_protspeps(protein_file, peptide_file, metadata_file)
#' 
#' SCsgp_species <- crumble_taxonomy(SC_subgroups, "species")
#' 
#' SCsgp_species_annot <- add_kegg(
#'   SCsgp_species, 
#'   kegg_db, 
#'   meta99_full_taxo, 
#'   metaproteome_origin, 
#'   protein_file, 
#'   peptide_file, 
#'   text_to_filter = "HUMAN"
#' )
#' 
#' }
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

add_kegg <-
  
function(spectral_count_object, annotation_db, taxonomic_db, 
         metaproteome_origin, protein_file, peptide_file, 
         text_to_filter = "HUMAN", taxonomic_levels_allowed = 1){

  genus <- Group <- ko <- number_ko <- Peptide <- phylum <- Protein <- species <- str_split <- Subgroup <- superkingdom <- NULL
  spectral_count_object <- spectral_count_object
  
  if (length(spectral_count_object) == 4){
    
    if (spectral_count_object[[4]] == "spectral_count_object"){
      
      pepsprots_origin <- spectral_count_object[[3]]
        
      # Check the format of protein_file
      proteins_info <- read.csv2(protein_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      if (is.null(proteins_info$SubGroup)){
        stop("Check proteins file format! A column named SubGroup should be present.")
      }
      else if (is.null(proteins_info$Protein)){
        stop("Check proteins file format! A column named Protein should be present.")
      }
      else if (is.null(proteins_info$Description)){
        stop("Check proteins file format! A column named Description should be present.")
      }
      else if (dim(proteins_info)[2] != 8){
        stop("Check proteins file format! Verify the number of columns.")
      }
      else{
        
        if(is.character(text_to_filter) == TRUE && nchar(text_to_filter) > 0){
          
          if("species.genus.family.order.class.phylum.superkingdom" %in% colnames(pepsprots_origin)){
            
            if("species" %in% colnames(pepsprots_origin) | 
               "genus" %in% colnames(pepsprots_origin) |
               "family" %in% colnames(pepsprots_origin) | 
               "order" %in% colnames(pepsprots_origin) | 
               "class" %in% colnames(pepsprots_origin) | 
               "phylum" %in% colnames(pepsprots_origin) | 
               "superkingdom" %in% colnames(pepsprots_origin)){
              
              annotation_db <- annotation_db
              
              if("gene_name" %in% colnames(annotation_db) & "ko" %in% colnames(annotation_db)){
                
                if(is.numeric(taxonomic_levels_allowed)){
                  if(round(taxonomic_levels_allowed, 0) >= 1){
                    
                    if((length(metaproteome_origin) == 6) && (metaproteome_origin[[6]] == "metaproteome_object")){
                      
                      # Check the format of peptide_file
                      peps_file <- read.csv2(peptide_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
                      if (is.null(peps_file$Subgroups)){
                        stop("Check peptides file format! A column named SubGroup should be present.")
                      }
                      else if (is.null(peps_file$Peptide)){
                        stop("Check peptides file format! A column named Peptide should be present.")
                      }
                      else if (is.null(peps_file$Group)){
                        stop("Check peptides file format! A column named Group should be present.")
                      }
                      else{
                        
                        proteins_info$Accession <- sapply(sapply(proteins_info$Description, str_split, pattern = " ", n = 2), "[[", 1)
                        taxonomic_db <- taxonomic_db
                        metaproteins_found <- subset(taxonomic_db, taxonomic_db[[1]] %in% proteins_info$Accession)
                        
                        # For functionnal annotation selection 
                        pepsprots <- merge(
                          proteins_info[, c("Group", "SubGroup", "Protein", "Description", "Accession")], 
                          metaproteins_found[, c(colnames(metaproteins_found)[1], "organism", "species.genus.family.order.class.phylum.superkingdom")],
                          by.x = "Accession",
                          by.y = colnames(metaproteins_found)[1],
                          all = TRUE
                        )
                        colnames(pepsprots)[3] <- "Subgroup"
                        pepsprots$organism <- as.character(pepsprots$organism)
                        pepsprots$species.genus.family.order.class.phylum.superkingdom <-
                          as.character(pepsprots$species.genus.family.order.class.phylum.superkingdom)
                        pepsprots[["species.genus.family.order.class.phylum.superkingdom"]][is.na(pepsprots[["species.genus.family.order.class.phylum.superkingdom"]])] <- "unassigned,unassigned,unassigned,unassigned,unassigned,unassigned,unassigned"
                        pepsprots[["organism"]][is.na(pepsprots[["organism"]])] <- "unassigned"
                        taxo_tab <- data.frame(do.call('rbind', strsplit(as.character(pepsprots$species.genus.family.order.class.phylum.superkingdom),',')))
                        colnames(taxo_tab) <- c("species", "genus", "family", "order", "class", "phylum", "superkingdom")
                        pepsprots <- cbind(pepsprots, taxo_tab)
                        
                        # Filter text from "Description" column
                        selected_elements <- pepsprots[grep(text_to_filter, pepsprots$Description), "Accession"]
                        pepsprots <- pepsprots[-which(pepsprots$Accession %in% selected_elements), ]
                        
                        # Taxonomic origin
                        if("species" %in% colnames(pepsprots_origin)){
                          taxa <- "species"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(species))), taxo_per_group = n_distinct(species))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(species))), taxo_per_group = n_distinct(species))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(species))), taxo_per_group = n_distinct(species))
                          }
                        }
                        else if("genus" %in% colnames(pepsprots_origin)){
                          taxa <- "genus"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(genus))), taxo_per_group = n_distinct(genus))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(genus))), taxo_per_group = n_distinct(genus))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(genus))), taxo_per_group = n_distinct(genus))
                          }
                        }
                        else if("family" %in% colnames(pepsprots_origin)){
                          taxa <- "family"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(family))), taxo_per_group = n_distinct(family))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(family))), taxo_per_group = n_distinct(family))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(family))), taxo_per_group = n_distinct(family))
                          }
                        }
                        else if("order" %in% colnames(pepsprots_origin)){
                          taxa <- "order"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(order))), taxo_per_group = n_distinct(order))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(order))), taxo_per_group = n_distinct(order))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(order))), taxo_per_group = n_distinct(order))
                          }
                        }
                        else if("class" %in% colnames(pepsprots_origin)){
                          taxa <- "class"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(class))), taxo_per_group = n_distinct(class))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(class))), taxo_per_group = n_distinct(class))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(class))), taxo_per_group = n_distinct(class))
                          }
                        }
                        else if("phylum" %in% colnames(pepsprots_origin)){
                          taxa <- "phylum"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(phylum))), taxo_per_group = n_distinct(phylum))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(phylum))), taxo_per_group = n_distinct(phylum))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(phylum))), taxo_per_group = n_distinct(phylum))
                          }
                        }
                        else{
                          taxa <- "superkingdom"
                          # Spectral origin
                          if(names(spectral_count_object[1]) == "SC_groups"){
                            sc_origin <- "Group"
                            taxo_assign <- pepsprots %>% group_by(Group) %>% summarize(taxa = names(which.max(table(superkingdom))), taxo_per_group = n_distinct(superkingdom))
                          }
                          else if(names(spectral_count_object[1]) == "SC_subgroups"){
                            sc_origin <- "Subgroup"
                            taxo_assign <- pepsprots %>% group_by(Subgroup) %>% summarize(taxa = names(which.max(table(superkingdom))), taxo_per_group = n_distinct(superkingdom))
                          }
                          else{
                            sc_origin <- "Peptide"
                            taxo_assign <- pepsprots %>% group_by(Protein) %>% summarize(taxa = names(which.max(table(superkingdom))), taxo_per_group = n_distinct(superkingdom))
                          }
                        }
                        taxo_assign <- as.data.frame(taxo_assign)
                        
                        # ---------- Quality Control plot
                        
                        # To keep graphics settings
                        old_par <- par(no.readonly = TRUE)
                        on.exit(suppressWarnings(par(old_par)))
                        
                        filename <- paste("qc", taxa, "per", sc_origin, sep = "_")
                        filename <- paste(filename, ".pdf", sep = "")
                        pdf(filename, width = 7.5, height = 6)
                        xlab <- paste(c("Number of ", sc_origin, "s"), collapse = "")
                        ylab <- paste(taxa, "per", sc_origin)
                        to_plot <- taxo_assign[order(taxo_assign[, 3]),]
                        plot(to_plot[, 3], pch = ".", cex = 2, ylab = ylab, xlab = xlab)
                        abline(h = 1, col = "red", lwd = 1)
                        dev.off()
                        print(paste("Quality control file", filename, "was generated"))
                        
                        # Selection of elements that contain only 
                        # 1 taxonomic level per spectral origin 
                        # (eg. 1 specie per subgroup)
                        names(taxo_assign)[2] <- taxa
                        total_groups <- dim(taxo_assign)[1]
                        selected_groups <- taxo_assign[taxo_assign$taxo_per_group <= taxonomic_levels_allowed, ]
                        nb_selected <- dim(selected_groups)[1]
                        taxa_percentage <- round(100 * (nb_selected / total_groups), 2)
                        print(paste("A total of", nb_selected, "out of", total_groups, "entities (", taxa_percentage, "% ) have a maximum of", taxonomic_levels_allowed, taxa, "per", sc_origin))
                        elements_to_keep <- as.vector(unique(selected_groups[[sc_origin]]))
                        new_pepsprots <- pepsprots[pepsprots[[sc_origin]] %in% elements_to_keep, ]
                        
                        # Integrate KEGG annotation
                        annotation_db$gene_name <- as.character(annotation_db$gene_name)
                        annotation_db$ko <- as.character(annotation_db$ko)
                        new_pepsprots$Accession <- as.character(new_pepsprots$Accession)
                        pepsprots_annotation <- merge(new_pepsprots, 
                                                      annotation_db[, c("gene_name", "ko")], 
                                                      by.x = "Accession", by.y = "gene_name", 
                                                      all.x = TRUE
                        )
                        pepsprots_annotation <- pepsprots_annotation[!is.na(pepsprots_annotation$ko), ]
                        nb_assigned_groups <- length(unique(as.factor(pepsprots_annotation[[sc_origin]])))
                        percent_ko_assigned <- round(100 * nb_assigned_groups / total_groups, 2)
                        print(paste("A total of", percent_ko_assigned, "% of", sc_origin, "s (", nb_assigned_groups, 
                                    "/", total_groups, ") were assigned to a KEGG Orthology term")
                        )
                        
                        # Kegg Orthology (KO) assignment, counts of KO per spectral
                        # groups or subgroups
                        print(paste0("Calculation of KO assignment per ", sc_origin, "s ...", collapse = ""))
                        
                        if(names(spectral_count_object[1]) == "SC_groups"){
                          ko_per_sc <- pepsprots_annotation %>% group_by(Group) %>% summarize(number_ko = n_distinct(ko), ko_uniques = paste(sort(unique(ko)), collapse = " "))
                          ko_per_prot <- pepsprots_annotation %>% group_by(Protein) %>% summarize(Group = names(which.max(table(Group))), ko_per_prot = n_distinct(ko), ko_uniques_prot = paste(sort(unique(ko)), collapse = " "))
                          ko_assign <- merge(ko_per_prot, ko_per_sc, by = "Group", all.x = TRUE)
                          prots_per_group <- pepsprots_annotation %>% group_by(Group) %>% summarize(prots_per_group = n_distinct(Protein))
                        }
                        else if(names(spectral_count_object[1]) == "SC_subgroups"){
                          ko_per_sc <- pepsprots_annotation %>% group_by(Subgroup) %>% summarize(number_ko = n_distinct(ko), ko_uniques = paste(sort(unique(ko)), collapse = " "))
                          ko_per_prot <- pepsprots_annotation %>% group_by(Protein) %>% summarize(Subgroup = names(which.max(table(Subgroup))), ko_per_prot = n_distinct(ko), ko_uniques_prot = paste(sort(unique(ko)), collapse = " "))
                          ko_assign <- merge(ko_per_prot, ko_per_sc, by = "Subgroup", all.x = TRUE)
                          prots_per_group <- pepsprots_annotation %>% group_by(Subgroup) %>% summarize(prots_per_group = n_distinct(Protein))
                        }
                        else{
                          ko_per_sc <- pepsprots_annotation %>% group_by(Protein) %>% summarize(number_ko = n_distinct(ko), ko_uniques = paste(sort(unique(ko)), collapse = " "))
                          ko_per_prot <- pepsprots_annotation %>% group_by(Protein) %>% summarize(Peptide = names(which.max(table(Peptide))), ko_per_prot = n_distinct(ko), ko_uniques_prot = paste(sort(unique(ko)), collapse = " "))
                          ko_assign <- merge(ko_per_prot, ko_per_sc, by = "Protein", all.x = TRUE)
                          prots_per_group <- pepsprots_annotation %>% group_by(Protein) %>% summarize(prots_per_group = n_distinct(Protein))
                        }
                        
                        # ---------- Quality Control plot
                        to_plot <- ko_assign
                        to_plot <- merge(to_plot, prots_per_group, by = sc_origin, all.x = TRUE)
                        filename <- paste("qc_KO_per", sc_origin, taxa, sep = "_")
                        filename <- paste(filename, ".pdf", sep = "")
                        xlab <- paste(c("Number of protein members per ", sc_origin), collapse = "")
                        ylab <- paste("Number of KO per", sc_origin)
                        p <- ggplot(to_plot, aes(x = prots_per_group, y = number_ko, alpha = 0.1)) + geom_jitter() + xlab(xlab) + ylab(ylab) + theme(legend.position = "none")
                        pdf(filename, width = 7.5, height = 6)
                        print(p)
                        dev.off()
                        print(paste("Quality control file", filename, "was generated"))
                        
                        # selection proteins and groups that fulfill the stated KO filters
                        ko_assign$selected_ko <- ifelse(ko_assign$number_ko == ko_assign$ko_per_prot, 1, 0)
                        ko_assign$union_ko <- ifelse(ko_assign$ko_uniques == ko_assign$ko_uniques_prot, 1, 0)
                        selected_elements <- ko_assign[ko_assign$selected_ko == 1 & ko_assign$union_ko == 1, ]
                        to_keep <- unique(selected_elements[[sc_origin]])
                        pepsprots_annotation <- pepsprots_annotation[pepsprots_annotation[[sc_origin]] %in% to_keep,]
                        pepsprots_annotation <- droplevels(pepsprots_annotation)
                        assigned_final <- length(unique(as.factor(pepsprots_annotation[[sc_origin]])))
                        percent_final <- round(100 * assigned_final / nb_assigned_groups, 2)
                        print(paste("From the", nb_assigned_groups, "elements assigned to a KO term,", assigned_final, 
                                    "(", percent_final, "% ) have redundant and consistent KO terms per", sc_origin)
                        )
                        
                        # Recalculate spectral counts after removal of proteins that do not meet the
                        # functional annotation requirements, method 2
                        removed_elements <- ko_assign[!ko_assign$selected_ko == 1 & !ko_assign$union_ko == 1, ]
                        removed_elements <- unique(removed_elements$Protein)
                        # Need to add peptides' id to remove the bad ones from the spectral data
                        pepsprots_annotation <- merge(
                          pepsprots_annotation, pepsprots_origin[, c("Protein", "Peptide")], 
                          by = "Protein", all.x = TRUE
                        )
                        removed_pepsprots <- pepsprots_annotation[pepsprots_annotation[["Protein"]] %in% removed_elements, ]
                        bad_peptides <- unique(removed_pepsprots$Peptide)
                        spectral_origin_data <- metaproteome_origin[[5]]
                        bad_sc_data <- spectral_origin_data[row.names(spectral_origin_data) %in% bad_peptides, ]
                        sc_to_remove <- merge(bad_sc_data, 
                                              unique(removed_pepsprots[, c(sc_origin, "Peptide", taxa)]),
                                              by.x = "row.names", 
                                              by.y = "Peptide"
                        )
                        sc_to_remove[[taxa]] <- as.character(sc_to_remove[[taxa]])
                        
                        #----- aggregated dataframe containing removed spectra
                        sc_data_removed <- apply(
                          sc_to_remove[, 2:(length(colnames(sc_to_remove))-2)], 
                          2, tapply, sc_to_remove[[taxa]], sum
                        )
                        sc_data_removed <- as.data.frame(sc_data_removed)
                        
                        #----- In case of any sample was removed during the analysis
                        metadata <- spectral_count_object[[2]]
                        msrun_to_keep <- metadata$SC_name
                        sc_data_removed <- sc_data_removed[, colnames(sc_data_removed) %in% msrun_to_keep]
                        zero_elements <- NULL
                        for (item in 1:dim(sc_data_removed)[1]){
                          all_zeros <- all(sc_data_removed[item,] == 0)
                          if (all_zeros == TRUE) {
                            zero_elements <- c(zero_elements, rownames(sc_data_removed[item,]))
                          }
                        }
                        sc_data_removed <- sc_data_removed[!rownames(sc_data_removed) %in% zero_elements,]
                        
                        #----- recalculation of spectral data, removal of spectra from taxonomic levels that do not met the filtering
                        pepsprots_annotation <- pepsprots_annotation[!is.na(pepsprots_annotation$Peptide), ]
                        taxo_to_keep <- as.vector(unique(pepsprots_annotation[[taxa]]))
                        sc_data <- spectral_count_object[[1]]
                        new_sc_data <- sc_data[rownames(sc_data) %in% taxo_to_keep, ]
                        
                        if(length(sc_data_removed) != 0){
                          decision_m1 <- stack(sc_data_removed[, 1:length(colnames(sc_data_removed))])
                          names(decision_m1)[1] <- "spectra_2_remove"
                          names(decision_m1)[2] <- "SC_name"
                          decision_m1 <- cbind.data.frame(
                            decision_m1, element = rep(rownames(sc_data_removed), length(colnames(sc_data_removed)))
                          )
                          names(decision_m1)[3] <- taxa
                          
                          decision_m2 <- stack(new_sc_data[, 1:length(colnames(new_sc_data))])
                          names(decision_m2)[1] <- "spectra"
                          names(decision_m2)[2] <- "SC_name"
                          decision_m2 <- cbind.data.frame(
                            decision_m2, element = rep(rownames(new_sc_data), length(colnames(new_sc_data)))
                          )
                          names(decision_m2)[3] <- taxa
                          
                          decision_matrix <- merge(decision_m2, decision_m1, by = c("SC_name", taxa), all.x = TRUE)
                          decision_matrix[is.na(decision_matrix)] <- 0
                          decision_matrix$new_spectra <- decision_matrix$spectra - decision_matrix$spectra_2_remove
                          
                          #----- aggregated dataframe containing updated spectra information
                          new_sc_data <- tapply(
                            decision_matrix[["new_spectra"]], 
                            list(decision_matrix[[taxa]], decision_matrix$SC_name), 
                            FUN=sum
                          )
                          new_sc_data <- as.data.frame(new_sc_data)
                        }
                        
                        #------- Reformat pepsprots dataframe
                        pepsprots_annotation <- merge(
                          pepsprots_annotation[, c("Accession", "Group", "Subgroup", "Protein", "Description", "Peptide", 
                                                   "organism", "species.genus.family.order.class.phylum.superkingdom", taxa, "ko")
                                               ],
                          pepsprots_origin[, c("Peptide", "Sequence", "Modifs", "MhTheo", "Charge")], 
                          by = "Peptide", all.x = TRUE
                        )
                        pepsprots_annotation <- pepsprots_annotation[, c("Accession", "Group", "Subgroup", "Protein", 
                                                                         "Description", "Peptide", "Sequence", "Modifs", 
                                                                         "MhTheo", "Charge", "organism", 
                                                                         "species.genus.family.order.class.phylum.superkingdom", 
                                                                         taxa, "ko")
                                                                     ]
                        pepsprots_annotation[[taxa]] <- as.character(pepsprots_annotation[[taxa]])
                        
                        #------------- Status message
                        print("========= Object information ============")
                        print(paste("Number of samples: ", length(unique(metadata$msrunfile)), sep = ""))
                        print(paste("Number of groups: ", length(unique(pepsprots_annotation$Group)), sep = ""))
                        print(paste("Number of subgroups: ", length(unique(pepsprots_annotation$Subgroup)), sep = ""))
                        print(paste("Number of peptides: ", length(unique(pepsprots_annotation$Peptide)), sep = ""))
                        print(paste("Number of ", taxa, " (taxonomic levels): ", dim(new_sc_data)[1], sep = ""))
                        print("=========================================")
                        output <- list(new_sc_data, metadata, pepsprots_annotation, "spectral_count_object")
                        names(output)[1] <- names(spectral_count_object[1])
                        names(output)[2] <- "metadata"
                        names(output)[3] <- "peptides_proteins"
                        names(output)[4] <- "type_object"
                        print("Spectral count object generated with the taxonomic elements having functional annotation")
                        return(output)
                      }
                    }
                    else{
                      stop("Invalid metaproteome_object, this object is generated from the function 'load_protspeps'")
                    }
                  }
                  else{
                    stop("The third argument must higher or equal to 1")
                  }
                }
                else{
                  stop("The third argument must be a numeric value")
                }
              }
              else{
                stop("Verify that the annotation DB has the variables: 'gene_name' and 'ko'")
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
          stop("The third argument must be a chain of characters")
        }
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
