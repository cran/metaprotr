#' export_robject
#' 
#' Exports one of the dataframes present in a "metaproteome_object"
#' or in "spectral_count_object". The export extensions can be RDATA or RDS.
#' 
#' @param entry_object A "metaproteome_object" or a "spectral_count_object" with
#'    similar format to that generated with the functions "load_protspeps" or 
#'    "getsc_specific", respectively.
#' 
#' @param data_exported Character indicating the type of data to be exported from 
#'    a "metaproteome_object" or a "spectral_count_object". The possible options 
#'    are: i) "proteins", ii) "peptides", iii) "pepProts", iv) "spectral", 
#'    v) "spectral_percent", vi) "metadata".
#' 
#' @param format_data Character indicating the file extension, this can be 
#'    either "RDATA" or "RDS". 
#'    
#' @param force Logic value set at FALSE by default in order to ask permission 
#'    to create object in the workstation of the user.
#' 
#' @return A file with the extension "RDATA" or "RDS" containing the information 
#'    from the selected dataframe from a "metaproteome_object" or a 
#'    "spectral_count_object".
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' export_robject(fecal_waters, "pepProts", "rdata")
#' 
#' data(species_fw)
#' export_robject(species_fw, "spectral", "rds")
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

export_robject <- 

function(entry_object, data_exported, format_data, force = FALSE){
  metaproteome_object <- entry_object
  
  if (!force && interactive()) {
    response <-
      select.list(c("yes", "no"), title = "Do you allow to create an R object?")
    
    if (response == "yes") {
      #-------- LOGIC TO EXPORT METAPROTEOME OBJECTS --------
      if (length(metaproteome_object) == 6) {
        
        if (metaproteome_object[[6]] == "metaproteome_object") {
          
          if (format_data == "rds" | format_data == "RDS") {
            
            if (data_exported == "proteins") {
              proteins_data <- metaproteome_object[[1]]
              save(proteins_data, file = "proteins_data.rds")
              print("'proteins_data.rds' file exported")
            }
            else if (data_exported == "peptides") {
              peptides_data <- metaproteome_object[[2]]
              save(peptides_data, file = "peptides_data.rds")
              print("'peptides_data.rds' file exported")
            }
            else if (data_exported == "pepProts") {
              peptides_proteins <- metaproteome_object[[3]]
              save(peptides_proteins, file = "peptides_proteins.rds")
              print("'peptides_proteins.rds' file exported")
            }
            else if (data_exported == "spectral") {
              spectral_counts <- metaproteome_object[[5]]
              save(spectral_counts, file = "spectral_counts.rds")
              print("'spectral_counts.rds' file exported")
            }
            else if (data_exported == "spectral_percent") {
              spectral_counts <- metaproteome_object[[5]]
              total_spectra <- as.data.frame(colSums(spectral_counts))
              total_spectra$SC_name <- rownames(total_spectra)
              names(total_spectra)[1] <- "total_spectra"
              sc_data_stacked <- stack(spectral_counts[, 1:length(colnames(spectral_counts))])
              names(sc_data_stacked)[1] <- "spectra"
              names(sc_data_stacked)[2] <- "SC_name"
              sc_data_stacked <- cbind.data.frame(sc_data_stacked, elements = rep(rownames(spectral_counts), length(colnames(spectral_counts))))
              sc_data_stacked[[2]] <- as.character(sc_data_stacked[[2]])
              sc_data_stacked[[3]] <- as.character(sc_data_stacked[[3]])
              sc_data_stacked <- merge(sc_data_stacked, total_spectra, by = "SC_name")
              sc_data_stacked$spectral_percent <- (sc_data_stacked$spectra / sc_data_stacked$total_spectra) * 100
              spectral_counts_percent <- tapply(
                sc_data_stacked$spectral_percent, 
                list(sc_data_stacked$elements, sc_data_stacked$SC_name), 
                FUN=sum
              )
              spectral_counts_percent <- as.data.frame(spectral_counts_percent)
              save(spectral_counts_percent, file = "spectral_counts_percent.rds")
              print("'spectral_counts_percent.rds' file exported")
            }
            else if (data_exported == "metadata") {
              metadata <- metaproteome_object[[4]]
              save(metadata, file = "metadata.rds")
              print("'metadata.rds' file exported")
            }
            else {
              stop(
                "Data exported can be ONE of the following options: i) 'proteins', ii) 'peptides', iii) 'pepProts', iv) 'spectral', iv) 'spectral_percent', vi) 'metadata'"
              )
            }
          }
          
          else if (format_data == "RDATA" | format_data == "rdata") {
            
            if (data_exported == "proteins") {
              proteins_data <- metaproteome_object[[1]]
              save(proteins_data, file = "proteins_data.RDATA")
              print("'proteins_data.RDATA' file exported")
            }
            else if (data_exported == "peptides") {
              peptides_data <- metaproteome_object[[2]]
              save(peptides_data, file = "peptides_data.RDATA")
              print("'peptides_data.RDATA' file exported")
            }
            else if (data_exported == "pepProts") {
              peptides_proteins <- metaproteome_object[[3]]
              save(peptides_proteins, file = "peptides_proteins.RDATA")
              print("'peptides_proteins.RDATA' file exported")
            }
            else if (data_exported == "spectral") {
              spectral_counts <- metaproteome_object[[5]]
              save(spectral_counts, file = "spectral_counts.RDATA")
              print("'spectral_counts.RDATA' file exported")
            }
            else if (data_exported == "spectral_percent") {
              spectral_counts <- metaproteome_object[[5]]
              total_spectra <- as.data.frame(colSums(spectral_counts))
              total_spectra$SC_name <- rownames(total_spectra)
              names(total_spectra)[1] <- "total_spectra"
              sc_data_stacked <- stack(spectral_counts[, 1:length(colnames(spectral_counts))])
              names(sc_data_stacked)[1] <- "spectra"
              names(sc_data_stacked)[2] <- "SC_name"
              sc_data_stacked <- cbind.data.frame(sc_data_stacked, elements = rep(rownames(spectral_counts), length(colnames(spectral_counts))))
              sc_data_stacked[[2]] <- as.character(sc_data_stacked[[2]])
              sc_data_stacked[[3]] <- as.character(sc_data_stacked[[3]])
              sc_data_stacked <- merge(sc_data_stacked, total_spectra, by = "SC_name")
              sc_data_stacked$spectral_percent <- (sc_data_stacked$spectra / sc_data_stacked$total_spectra) * 100
              spectral_counts_percent <- tapply(
                sc_data_stacked$spectral_percent, 
                list(sc_data_stacked$elements, sc_data_stacked$SC_name), 
                FUN=sum
              )
              spectral_counts_percent <- as.data.frame(spectral_counts_percent)
              save(spectral_counts_percent, file = "spectral_counts_percent.RDATA")
              print("'spectral_counts_percent.RDATA' file exported")
            }
            else if (data_exported == "metadata") {
              metadata <- metaproteome_object[[4]]
              save(metadata, file = "metadata.RDATA")
              print("'metadata.RDATA' file exported")
            }
            else {
              stop(
                "Data exported can be ONE of the following options: i) 'proteins', ii) 'peptides', iii) 'pepProts', iv) 'spectral', iv) 'spectral_percent', vi) 'metadata'"
              )
            }
          }
          
          else{
            stop("Only RDATA or rds formats allowed, select one of them")
          }
          
        }
        
        else{
          stop("Invalid metaproteome object")
        }
      }
      
      #-------- LOGIC TO EXPORT SPECTRAL OBJECTS --------
      else if (length(metaproteome_object) == 4) {
        
        if (metaproteome_object[[4]] == "spectral_count_object") {
          
          if (format_data == "rds" | format_data == "RDS") {
            if (data_exported == "spectral") {
              spectral_specific <- metaproteome_object[[1]]
              saveRDS(spectral_specific, file = "spectral_counts.rds")
              print("'spectral_specific.rds' file exported")
            }
            else if (data_exported == "spectral_percent") {
              spectral_counts <- metaproteome_object[[1]]
              total_spectra <- as.data.frame(colSums(spectral_counts))
              total_spectra$SC_name <- rownames(total_spectra)
              names(total_spectra)[1] <- "total_spectra"
              sc_data_stacked <- stack(spectral_counts[, 1:length(colnames(spectral_counts))])
              names(sc_data_stacked)[1] <- "spectra"
              names(sc_data_stacked)[2] <- "SC_name"
              sc_data_stacked <- cbind.data.frame(sc_data_stacked, elements = rep(rownames(spectral_counts), length(colnames(spectral_counts))))
              sc_data_stacked[[2]] <- as.character(sc_data_stacked[[2]])
              sc_data_stacked[[3]] <- as.character(sc_data_stacked[[3]])
              sc_data_stacked <- merge(sc_data_stacked, total_spectra, by = "SC_name")
              sc_data_stacked$spectral_percent <- (sc_data_stacked$spectra / sc_data_stacked$total_spectra) * 100
              spectral_counts_percent <- tapply(
                sc_data_stacked$spectral_percent, 
                list(sc_data_stacked$elements, sc_data_stacked$SC_name), 
                FUN=sum
              )
              spectral_counts_percent <- as.data.frame(spectral_counts_percent)
              save(spectral_counts_percent, file = "spectral_counts_percent.rds")
              print("'spectral_counts_percent.rds' file exported")
            }
            else if (data_exported == "metadata") {
              metadata <- metaproteome_object[[2]]
              saveRDS(metadata, file = "metadata.rds")
              print("'metadata.rds' file exported")
            }
            else if (data_exported == "pepProts") {
              peptides_proteins <- metaproteome_object[[3]]
              save(peptides_proteins, file = "peptides_proteins.rds")
              print("'peptides_proteins.rds' file exported")
            }
            
            else {
              stop(
                "Data exported can be ONE of the following options: i)'spectral', ii) 'spectral_percent', iii) 'metadata', iv) 'pepProts'"
              )
            }
          }
          
          else if (format_data == "RDATA" | format_data == "rdata") {
            
            if (data_exported == "spectral") {
              spectral_specific <- metaproteome_object[[1]]
              save(spectral_specific, file = "spectral_counts.RDATA")
              print("'spectral_specific.RDATA' file exported")
            }
            else if (data_exported == "spectral_percent") {
              spectral_counts <- metaproteome_object[[1]]
              total_spectra <- as.data.frame(colSums(spectral_counts))
              total_spectra$SC_name <- rownames(total_spectra)
              names(total_spectra)[1] <- "total_spectra"
              sc_data_stacked <- stack(spectral_counts[, 1:length(colnames(spectral_counts))])
              names(sc_data_stacked)[1] <- "spectra"
              names(sc_data_stacked)[2] <- "SC_name"
              sc_data_stacked <- cbind.data.frame(sc_data_stacked, elements = rep(rownames(spectral_counts), length(colnames(spectral_counts))))
              sc_data_stacked[[2]] <- as.character(sc_data_stacked[[2]])
              sc_data_stacked[[3]] <- as.character(sc_data_stacked[[3]])
              sc_data_stacked <- merge(sc_data_stacked, total_spectra, by = "SC_name")
              sc_data_stacked$spectral_percent <- (sc_data_stacked$spectra / sc_data_stacked$total_spectra) * 100
              spectral_counts_percent <- tapply(
                sc_data_stacked$spectral_percent, 
                list(sc_data_stacked$elements, sc_data_stacked$SC_name), 
                FUN=sum
              )
              spectral_counts_percent <- as.data.frame(spectral_counts_percent)
              save(spectral_counts_percent, file = "spectral_counts_percent.RDATA")
              print("'spectral_counts_percent.RDATA' file exported")
            }
            else if (data_exported == "metadata") {
              metadata <- metaproteome_object[[2]]
              save(metadata, file = "metadata.RDATA")
              print("'metadata.RDATA' file exported")
              
            }
            else if (data_exported == "pepProts") {
              peptides_proteins <- metaproteome_object[[3]]
              save(peptides_proteins, file = "peptides_proteins.RDATA")
              print("'peptides_proteins.RDATA' file exported")
              
            }
            else {
              stop(
                "Data exported can be ONE of the following options: i)'spectral', ii) 'spectral_percent', iii) 'metadata', 'iv) pepProts'"
              )
            }
          }
          
          else{
            stop("Only RDATA or rds formats allowed, select one of them")
          }
        }
        
        else{
          stop("Invalid spectral count object")
        }
      }
      
      else{
        stop("Invalid object")
      }
      
    }
    else{
      stop("No file was created")
    }
  }      
}
