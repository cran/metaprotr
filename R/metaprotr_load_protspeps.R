#' load_protspeps
#' 
#' Loads three files: i) peptides abundances expressed as 
#' spectral counts, ii) proteins information, and iii) metadata of the
#' mass spectrometry samples. Combines the three files into a 
#' "metaproteome_object", a list containing these dataframes.
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
#' @param metadata_file Character indicating the location of a csv file 
#'    containing the samples information. The following columns names 
#'    MUST be present: "SC_name" (sample ids assigned by the user), "msrunfile" 
#'    (name of samples as indicated in mass spectrometry files and in the columns
#'    of peptide_file) and "SampleID" (codes indicating the experimental group).
#'    Additional columns containing complementary information can be added by 
#'    the user (ex. replicates, order of injection, etc.). 
#'    Separation between columns should be indicated by tabulation. 
#'    For more details regarding data input check 
#'    \href{https://forgemia.inra.fr/pappso/metaprotr#data-inputs}{format examples}.
#'
#' @return A "metaproteome_object", which is a list of six elements containing: 
#'    1) dataframe of the protein identifiers, 
#'    2) dataframe of the peptide identifiers, 
#'    3) dataframe containing the information of peptides with their associated proteins, 
#'    4) dataframe of metadata containing the experiment information, 
#'    5) dataframe of spectral counts per peptide on each sample, 
#'    6) character indicating the type of object generated. 
#'
#' @examples
#' \dontrun{
#' 
#' protein_file <- "location/peptides_abundances.csv"
#' peptide_file <- "location/proteins_list.csv"
#' metadata <- "location/metadata.csv"
#' metaproteome <- load_protspeps(protein_file, peptide_file, metadata_file)
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

load_protspeps <- 

function(protein_file, peptide_file, metadata_file){
  
  str_split <- NULL
  proteins_info <- read.csv2(protein_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  peps_file <- read.csv2(peptide_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  metadata <- read.csv2(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Conditionals of the correct file format of protein_file
  # Column 'SubGroup' is found on proteins list file
  if (is.null(proteins_info$SubGroup)) {
    stop("Check proteins file format! A column named SubGroup should be present.")
  }
  
  else if (is.null(proteins_info$Protein)) {
    stop("Check proteins file format! A column named Protein should be present.")
  }
  
  else if (is.null(proteins_info$Description)) {
    stop("Check proteins file format! A column named Description should be present.")
  }
  
  else if (dim(proteins_info)[2] != 8) {
    stop("Check proteins file format! Verify the number of columns.")
  }
  
  # Conditionals to check the correct file format of peptide_file
  # Column "Subgroups" is found in petides file
  else if (is.null(peps_file$Subgroups)) {
    stop("Check peptides file format! A column named SubGroup should be present.")
  }
  
  else if (is.null(peps_file$Peptide)) {
    stop("Check peptides file format! A column named Peptide should be present.")
  }
  
  else if (is.null(peps_file$Group)) {
    stop("Check peptides file format! A column named Group should be present.")
  }
  
  # Conditionals of the correct file format of metadata_file
  else if (is.null(metadata$msrunfile)) {
    stop("Check metadata file format! A column named msrunfile should be present. The columns must be separated by tabulation.")
  }
  
  else if (is.null(metadata$SampleID)) {
    stop("Check metadata file format! A column named SampleID should be present. The columns must be separated by tabulation.")
  }
  
  else if (is.null(metadata$SC_name)) {
    stop("Check metadata file format! A column named SC_name should be present. The columns must be separated by tabulation.")
  }
  
  else{
    
    # ------------------ LOAD PROTEINS FILE ------------------
    print("Loading proteins...")
    dimension_origin <- dim(proteins_info)
    length_origin <- length(unique(proteins_info$SubGroup))
    
    # keep proteins finishing with ".a1" 
    proteins_info <- proteins_info[grepl(".a1$", proteins_info$Protein), ]
    
    # Assign "accession" identifier from the first string of characters in the "Description" column
    proteins_info$Accession <-
      sapply(
        sapply(proteins_info$Description, str_split, pattern = " ", n = 2), "[[", 1
      )
    dimension_final <- dim(proteins_info)
    
    # print(paste(c("Proteins file with", dimension_final[1], "proteins in", dimension_final[2], "columns"), collapse = " ") )
    # ------------------ LOAD PEPTIDES FILE ------------------
    
    print("Loading peptides...")
    # Total number of peptides
    total_peps <- dim(peps_file)[1]
    all_peps_info <- peps_file[, c(1:7)]
    
    if (dim(all_peps_info)[2] != 7){
      stop("Check the peptides file format! Should have 7 columns with peptides metadata.")
    }
    
    else {
      all_peps_list <- strsplit(as.character(all_peps_info[, "Subgroups"]), " ")
      
      # Selection of specific peptide
      specific_peps <- lapply(all_peps_list, length) == 1
      pepdata <- all_peps_info[specific_peps, c("Group", "Subgroups", "Peptide")]
      count_specific_peps <- dim(pepdata)[1]
      
      # Selection of non-specific peptides
      shared_peps <- all_peps_info[!specific_peps, ]
      count_shared_peps <- dim(shared_peps)[1]
      
      # List of dataframes containing the groups and subgroups per peptide
      shared_peps_list <- strsplit(as.character(shared_peps[, "Subgroups"]), " ")
      nested_shared_peps <-
        lapply(shared_peps_list, 
               function(x) {data.frame(Group = gsub("\\.[a-z]+[0-9]+$", "", x), Subgroups = x)})
      names(nested_shared_peps) <- shared_peps$Peptide
      
      # Add peptide ID per dataframe
      for (i in 1:length(nested_shared_peps)) {
        # if(i%%5000 == 0){
        #   print(i)
        # }
        nested_shared_peps[[i]]$Peptide <-
          rep(names(nested_shared_peps)[[i]], nrow(nested_shared_peps[[i]]))
      }
      print("This step takes some time...")
      nested_shared_peps[[length(nested_shared_peps) + 1]] <- pepdata
      pepdata <- do.call("rbind", nested_shared_peps)
      print(paste(
        c("Data with", count_specific_peps, "specific peptides and",
          count_shared_peps, "non-specific peptides"), collapse = " ")
      )
      
      # ------------------ COMBINE PROTEINS AND PEPTIDES ------------------
      pep_prot <- merge(proteins_info, pepdata, by.x = "SubGroup", by.y = "Subgroups")
      pep_prot <-pep_prot[, c("Group.x", 
                              "SubGroup", 
                              "Protein", 
                              "Accession", 
                              "Description", 
                              "Peptide")
                          ]
      colnames(pep_prot)[1] <- "Group"
      
      # merge [link subgroup - name of protein .a1 - peptide] and infos on peptides
      pep_prot <- merge(pep_prot, all_peps_info[, c("Peptide", "Sequence", "Modifs", "MhTheo", "Charge")], by = "Peptide")
      # Reorder
      pep_prot <- pep_prot[, c(
        "Group", "SubGroup", "Protein", "Accession", "Description", "Peptide", 
        "Sequence", "Modifs", "MhTheo", "Charge")
        ]
      colnames(pep_prot)[2] <- "Subgroup"
      print("Metadata file loaded")
      
      # ------------------ SPECTRAL COUNTING OF PEPTIDES ------------------
      
      rownames(peps_file) <- peps_file$Peptide
      SCpep_matrix <- peps_file[, 8:(dim(peps_file)[2])]
      
      # Checks if that msrunfile starts with a numerical value
      if (startsWith(names(SCpep_matrix)[1], "X") == TRUE) {
        colnames(SCpep_matrix) <- gsub("^X", "" , colnames(SCpep_matrix))
      }
      
      # Content of "msrunfile" matches that of "peptide spectral counting"
      if (all(colnames(as.character(SCpep_matrix)) %in% metadata$msrunfile) == FALSE) {
        stop(
          "Verify that 'msrunfiles' content is the same in metadata and peptides_counting files"
        )
      }
      else{
        # SCpep_matrix <- SCpep_matrix[, metadata$msrunfile]
        SCpep <- SCpep_matrix
        colnames(SCpep) <- metadata$SC_name
        print("Spectral counting of peptides added")
        
        # Status message
        print("========= Object information ============")
        print(paste("Number of samples: ", length(unique(metadata$msrunfile)), sep = ""))
        print(paste("Number of groups: ", length(unique(pep_prot$Group)), sep = ""))
        print(paste("Number of subgroups or metaproteins: ", length(unique(pep_prot$Subgroup)), sep = ""))
        print(paste("Number of peptides: ", length(unique(pep_prot$Peptide)), sep = ""))
        print("=========================================")
        
        # OUTPUT DATA 
        output <-
          list(proteins_info,
               pepdata,
               pep_prot,
               metadata,
               SCpep,
               "metaproteome_object")
        names(output)[1] <- "proteins_data"
        names(output)[2] <- "peptides_data"
        names(output)[3] <- "peptides_proteins"
        names(output)[4] <- "metadata"
        names(output)[5] <- "spectral_counting_peptides"
        names(output)[6] <- "type_object"
        print("metaproteome object generated")
        return(output)
      }
    }
  }
}