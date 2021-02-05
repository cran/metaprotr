#' plot_stackedtaxo
#'
#' Generates stacked barplots of the spectral counts distributions among the different
#' taxonomic entities ("species", "genus", "family", "order", "class", "phylum" 
#' or "superkingdom") within the samples or conditions of a "spectral_count_object" 
#' with taxonomy. If the provided conditions have several replicates 
#' the mean value is taken into account.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with spectral counts abundance organized by taxonomy 
#'    (species, genus, family, order, class, phylum or superkingdom). This 
#'    object is generated with the function "crumble_taxonomy".
#'
#' @param target_variable Character indicating the name of one column from 
#'    metadata. The stacked barplots will be ordered by the levels of this variable.
#'
#' @param bars_data Character indicating the type of labels to be displayed in 
#'    the stacked bars. The possible options are "percent" or "numbers".
#' 
#' @param filter_percent Optional numeric value between 0 and 99 that 
#'    sets the minimal percentage of spectral counts at which the taxonomic 
#'    elements will be displayed. The elements whose values are lower than 
#'    this number will be gathered and displayed as "others". 
#'    The default value is set at 1.
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a pdf in the workstation of the user.
#'
#' @return Barplots (pdf) of the taxonomic distribution of the samples present 
#'    in a "spectral_count_object" with taxonomic levels. 
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(species_fw)
#' 
#' plot_stackedtaxo(species_fw, 'SampleID', 'percent', 2)
#' 
#' plot_stackedtaxo(species_fw, 'SC_name', 'numbers')
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

plot_stackedtaxo <-

function(spectral_count_object, target_variable, 
         bars_data, filter_percent = 1, force = FALSE){
  
  percent_spectra <- spectra <- taxonomy <- NULL
  spectral_count_object <- spectral_count_object
  
  if(!force && interactive()){
    response <- select.list(
      c("yes", "no"),
      title = paste("Do you allow to create a plot with taxonomic distributions per", target_variable, "?", sep = " ")
    )
    if(response == "yes"){
      if (length(spectral_count_object) == 4) {
        if(spectral_count_object[[4]] == "spectral_count_object") {
          
          pepsprots <- spectral_count_object[[3]]
          
          if("species.genus.family.order.class.phylum.superkingdom" %in% colnames(pepsprots)){
            if("species" %in% colnames(pepsprots) | "genus" %in% colnames(pepsprots) |
               "family" %in% colnames(pepsprots) | "order" %in% colnames(pepsprots) | 
               "class" %in% colnames(pepsprots) | "phylum" %in% colnames(pepsprots) | 
               "superkingdom" %in% colnames(pepsprots)){
              
              if("species" %in% colnames(pepsprots)){
                taxa <- "species"
              }
              else if("genus" %in% colnames(pepsprots)){
                taxa <- "genus"
              }
              else if("family" %in% colnames(pepsprots)){
                taxa <- "family"
              }
              else if("order" %in% colnames(pepsprots)){
                taxa <- "order"
              }
              else if("class" %in% colnames(pepsprots)){
                taxa <- "class"
              }
              else if("phylum" %in% colnames(pepsprots)){
                taxa <- "phylum"
              }
              else{
                taxa <- "superkingdom"
              }
              metadata <- spectral_count_object[[2]]
              
              if(target_variable %in% colnames(metadata)){
                if(bars_data == "percent" | bars_data == "numbers"){
                  if(filter_percent < 0 | filter_percent > 99){
                    stop("The fourth argument must be a number between 0 and 99")
                  }
                  else{
                    sc_data <- spectral_count_object[[1]]
                    colnames(sc_data) <- metadata[[target_variable]]
                    sc_data <- sapply(split.default(sc_data, names(sc_data)), rowMeans)
                    sc_data <- as.data.frame(sc_data)
                    
                    # To display the top x% elements
                    total_spectra <- as.data.frame(colSums(sc_data))
                    total_spectra[[target_variable]] <- rownames(total_spectra)
                    names(total_spectra)[1] <- "total_spectra"
                    sc_data_stacked <- stack(sc_data[, 1:length(colnames(sc_data))])
                    names(sc_data_stacked)[1] <- "spectra"
                    names(sc_data_stacked)[2] <- target_variable
                    sc_data_stacked <- cbind.data.frame(
                      sc_data_stacked, 
                      taxo_levels = rep(rownames(sc_data), length(colnames(sc_data))))
                    names(sc_data_stacked)[3] <- taxa
                    sc_data_stacked[[2]] <- as.character(sc_data_stacked[[2]])
                    sc_data_stacked[[3]] <- as.character(sc_data_stacked[[3]])
                    sc_data_stacked <- merge(sc_data_stacked, total_spectra, by=target_variable)
                    sc_data_stacked$percent_spectra <- round(sc_data_stacked$spectra / sc_data_stacked$total_spectra, 5) * 100
                    
                    # Filter provided by the user
                    sc_data_stacked$taxonomy <- ifelse(sc_data_stacked$percent_spectra >= filter_percent, sc_data_stacked[[taxa]], "Others")
                    sc_data_stacked <- sc_data_stacked %>% group_by(
                      sc_data_stacked[[target_variable]], taxonomy) %>% summarize(
                        spectra = round(sum(spectra), 0), 
                        percent_spectra = round(sum(percent_spectra), 3))
                    names(sc_data_stacked)[1] <- target_variable
                    sc_data_stacked <- as.data.frame(sc_data_stacked)
                    
                    # Plot logic
                    # To keep graphics settings
                    old_par <- par(no.readonly = TRUE)
                    on.exit(suppressWarnings(par(old_par)))
                    
                    if(bars_data == "numbers"){
                      variable_x <- sc_data_stacked[[target_variable]]
                      ylab <- paste("Number of spectra per", taxa, sep = " ")
                      p <- ggplot(
                        sc_data_stacked, aes(fill = taxonomy, y = spectra, x = variable_x)) + 
                        geom_bar(position="stack", stat="identity") +
                        labs(x = "", y = ylab) +
                        scale_fill_viridis_d(direction = -1) +
                        theme_classic(base_size = 11) +
                        theme(legend.position = "bottom")
                    }
                    else{
                      variable_x <- sc_data_stacked[[target_variable]]
                      ylab <- paste("Percentage of spectra per", taxa, sep = " ")
                      p <- ggplot(
                        sc_data_stacked, aes(fill = taxonomy, y = percent_spectra, x = variable_x)) + 
                        geom_bar(position="stack", stat="identity") +
                        labs(x = "", y = ylab) +
                        scale_fill_viridis_d(direction = -1) +
                        theme_classic(base_size = 11) +
                        theme(legend.position = "bottom")
                    }
                    filename <- paste("Bars_taxonomy", taxa, target_variable, bars_data, filter_percent, sep = "_")
                    filename <- paste(filename, ".pdf", sep = "")
                    pdf(filename, width = 16, height = 14)
                    print(p)
                    dev.off()
                    print(paste("The figure", filename, "was generated", sep = " "))
                  }
                }
                else{
                  stop("The third argument must be: 'percent' or 'numbers'")
                }
              }
              else{
                var_options <- colnames(metadata)
                stop(paste(
                  c("The second argument must be one of the following options", 
                    var_options), 
                  collapse = ", "))
              }
            }
            else{
              stop("Invalid spectral_count_object, spectral counts must be reorganized with the function crumble_taxonomy")
            }
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