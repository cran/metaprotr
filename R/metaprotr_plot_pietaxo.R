#' plot_pietaxo
#'
#' Generates a pie chart with taxonomic distribution of one selected
#' sample or condition. If the provided condition has several replicates 
#' the mean value is taken into account.
#' 
#' @param spectral_count_object list defined as "spectral_count_object" 
#'    containing dataframes with spectral counts abundance of the 
#'    samples organized by taxonomy (species, genus, family, 
#'    order, class, phylum or superkingdom). This object is generated 
#'    with the function "crumble_taxonomy".
#'
#' @param target_variable Character indicating the name of one column 
#'    from metadata. This column must contain the identifiers of the sample 
#'    or condition to be plotted.
#'
#' @param sampling Character indicating the name of sample or condition
#'    to be plotted. This character must be present in the "target_variable". 
#'  
#' @param filter_percent Optional numeric value between 0 and 99 that 
#'    sets the minimal percentage of spectral counts at which the taxonomic 
#'    elements will be displayed. The elements whose values are lower than 
#'    this number will be gathered and displayed as "others". 
#'    The default value is set at 1.
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create files in the workstation of the user.
#'
#' @return A pie chart (pdf) and a csv file with the taxonomic 
#'    distribution of one sample or one condition. In the csv file, 
#'    all the elements have their: a) spectral counts, b) percentage of 
#'    spectral count in the sample, c) the taxonomic name and d)
#'    the taxonomic elements that were assiged as 'others' in function 
#'    of the filter provided. 
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(species_fw)
#' 
#' plot_pietaxo(species_fw, "Methods", "S")
#' 
#' plot_pietaxo(species_fw, "SC_name", "Q1")
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

plot_pietaxo <-

function(spectral_count_object, target_variable, 
         sampling, filter_percent = 1, force = FALSE){
  
  percent_spectra <- spectra <- taxonomy <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <- select.list(
      c("yes", "no"),
      title = paste("Do you allow to create a pie chart and csv file with taxonomic distributions per", target_variable, "?", sep = " ")
    )
    if (response == "yes") {
      if (length(spectral_count_object) == 4) {
        if(spectral_count_object[[4]] == "spectral_count_object") {
          
          pepsprots <- spectral_count_object[[3]]
          
          if("species.genus.family.order.class.phylum.superkingdom" %in% colnames(pepsprots)){
            if("species" %in% colnames(pepsprots) | "genus" %in% colnames(pepsprots) |
               "family" %in% colnames(pepsprots) | "order" %in% colnames(pepsprots) | 
               "class" %in% colnames(pepsprots) | "phylum" %in% colnames(pepsprots) | 
               "superkingdom" %in% colnames(pepsprots)){
              # For aesthetics in plot
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
                options_var <- metadata[[target_variable]]
                if(sampling %in% options_var){
                  if(filter_percent < 0 | filter_percent > 99){
                    stop("The fourth argument must be a number between 0 and 99")
                  }
                  else{
                    # data wrangling
                    sc_data <- spectral_count_object[[1]]
                    colnames(sc_data) <- metadata[[target_variable]]
                    sc_data <- round(sapply(split.default(sc_data, names(sc_data)), rowMeans), 2)
                    sc_data <- as.data.frame(sc_data)
                    sc_data <- sc_data[which(sc_data[[sampling]] > 0), sampling, drop = FALSE]
                    total_spectra <- sum(sc_data[[sampling]])
                    sc_data$percent_spectra <- 100 * (sc_data[[sampling]] / total_spectra)
                    sc_data[[taxa]] <- rownames(sc_data)
                    names(sc_data)[1] <- "spectra"
                    
                    # To display the top x% elements
                    sc_data$taxonomy <- ifelse(sc_data$percent_spectra >= filter_percent, sc_data[[taxa]], "Others")
                    sc_data$taxonomy <- as.factor(sc_data$taxonomy)
                    others_count <- length(which(sc_data$taxonomy == 'Others'))
                    sc_data_summary <- as.data.frame(tapply(sc_data$spectra, sc_data$taxonomy, sum))
                    names(sc_data_summary)[1] <- "spectra"
                    sc_data_summary$taxonomy <- rownames(sc_data_summary)
                    sc_data_summary$taxonomy <- ifelse(
                      sc_data_summary$taxonomy == "Others", 
                      paste0("Others [", others_count, " entities]"),
                      sc_data_summary$taxonomy
                    )
                    sc_data_summary$taxonomy <- as.factor(sc_data_summary$taxonomy)
                    sc_data_summary$percent_spectra <- round(100 * (sc_data_summary$spectra / total_spectra), 2)
                    sc_data_summary$spectra <- round(sc_data_summary$spectra, 0)
                    
                    # To handle lots of entities
                    if(filter_percent < 1){
                      sc_data_summary$label_plot <- paste0(
                        sc_data_summary$spectra, 
                        " (", 
                        round(sc_data_summary$percent_spectra, 1), 
                        " %)\n", 
                        sc_data_summary$taxonomy
                      )
                    }
                    else{
                      sc_data_summary$label_plot <- ifelse(
                        as.numeric(sc_data_summary$percent_spectra) >= 1, 
                        paste0(
                          sc_data_summary$spectra, 
                          " (", 
                          round(sc_data_summary$percent_spectra, 1), 
                          " %)\n", 
                          sc_data_summary$taxonomy
                        ), 
                        " "
                      )
                    }
                    
                    # Plot logic
                    # To keep graphics settings
                    old_par <- par(no.readonly = TRUE)
                    on.exit(suppressWarnings(par(old_par)))
                    
                    if(taxa == "genus"){
                      taxa <- "genera"
                    }
                    main <- paste("Distribution of spectral counts among", taxa, "in", sampling, sep = " ")
                    p <- ggplot(sc_data_summary, aes(x = "", y = spectra, fill = taxonomy)) +
                      geom_bar(stat="identity", width = 1, color = "white") +
                      scale_fill_viridis_d(direction = -1) +
                      geom_text_repel(
                        aes(x = 2, label = label_plot),
                        position = position_stack(vjust = 0.01), 
                        color = "black") +
                      labs(x = "", title = main) +
                      coord_polar("y", start = 0) +
                      theme_void(base_size = 10) +
                      theme(legend.position = "bottom")
                    filename <- paste("Pie", taxa, target_variable, sampling, sep = "_")
                    filename_csv <- paste(filename, ".csv", sep = "")
                    filename <- paste(filename, ".pdf", sep = "")
                    pdf(filename, width = 10, height = 19)
                    print(p)
                    dev.off()
                    write.table(sc_data, filename_csv, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
                    print(paste(filename, "and", filename_csv, "were generated", sep = " "))
                  }
                }
                else{
                  vars <- unique(options_var)
                  stop(paste(c("The third argument must be one these options", vars), collapse = ", "))
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