#' plot_dendocluster
#'
#' Draws a dendogram where samples are clustered based on the number of elements 
#' present on each sample from a "spectral_count_object". This graph is 
#' constructed based on Spearman correlations transformed into distances 
#' and plotted with the logic of the package 
#' \href{https://CRAN.R-project.org/package=dendextend}{dendextend}.
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts 
#'    from peptides, subgroups, groups or taxonomic levels. The format 
#'    of this object is similar to that generated from the functions
#'    "getsc_specific" and "crumble_taxonomy".
#' 
#' @param target_variable Character indicating the name of one column 
#'    from metadata. The different levels in this column will be represented 
#'    as different colors in the final dendogram.
#' 
#' @param file_title Character indicating the name of the generated file.
#'  
#' @param hclust_method Character indicating the agglomeration method to be 
#'    used for the hierarchical clustering. The possible methods are discribed 
#'    on \href{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust}{hclust}.
#'    The default method is "ward.D".
#  
#' @param correlation_method Character indicating the correlation coeficient to 
#'    be computed. The possible options are discribed in the function 
#'    \href{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor}{cor}. 
#'    The default value is "spearman".
#' 
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a pdf file in the workstation of the user. 
#'
#' @return A dendogram plot (pdf) indicating the number of elements per
#'    sample.
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' str(fecal_waters$metadata)
#' 
#' plot_dendocluster(fecal_waters, "Condition", "title_dendogram")
#' 
#' plot_dendocluster(fecal_waters, "Condition", "title_dendogram", 
#'    hclust_method = "mcquitty")
#' 
#' plot_dendocluster(fecal_waters, "Condition", "title_dendogram_groups", 
#'    correlation_method = "pearson")
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

plot_dendocluster <- 

function(spectral_count_object, target_variable, file_title, 
         hclust_method = "ward.D", correlation_method = "spearman", force = FALSE){
  
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    
    response <- select.list(c("yes", "no"), title = "Do you allow to create a pdf file with a dendogram?")
    
    if (response == "yes") {
      if ((length(spectral_count_object) == 4) &
          (spectral_count_object[[4]] == "spectral_count_object")) {
        
        metadata <- spectral_count_object[[2]]
        
        if(target_variable %in% colnames(metadata)){
          
          elements_meta <- as.character(metadata[[target_variable]])
          # list of colors
          colors_condition = NULL
          unique_conditions <- unique(elements_meta)
          print("Indicate colors by name (ex. red) or hexadecimal code (ex. #56338C)")
          for(i in 1:length(unique_conditions)){
            color <- readline(prompt = paste0("Enter the color for ", unique_conditions[i], ": " ))
            colors_condition <- rbind.data.frame(
              colors_condition, 
              cbind.data.frame(condition = as.character(unique_conditions[i]), color = as.character(color))
            )
          }
          # assign ordered colors
          colors_vector = NULL
          for(i in 1:length(elements_meta)){
            for(j in 1:dim(colors_condition)[1]){
              if(elements_meta[i] == as.character(colors_condition[j, 1])){
                colors_vector <- append(colors_vector, as.character(colors_condition[j, 2]))
              }
            }
          }
          ncols <- length(colnames(spectral_count_object[[1]]))
          ncolors <- length(colors_vector)
          
          if (ncols == ncolors) {
            if (is.character(file_title) == TRUE && nchar(file_title) > 1){
              SC_sampple_id <- spectral_count_object[[1]]
              colnames(SC_sampple_id) <- metadata$SampleID
              
              # Add the number of elements per sample
              max_length_name <- max(nchar(colnames(SC_sampple_id)))
              for (sample in 1:c(length(colnames(SC_sampple_id)))) {
                sample_name <- names(SC_sampple_id[sample])
                length_name <- nchar(sample_name)
                length_space <- max_length_name - length_name + 7
                length_space <- paste(replicate(length_space, " "), collapse = "")
                sample_elements <- nrow(SC_sampple_id[SC_sampple_id[sample] > 0, ])
                new_name <- paste(sample_name, length_space, sample_elements)
                names(SC_sampple_id)[sample] <- new_name
              }
              
              # Spearman correlation
              SC_sampple_id_corr <- 1 - cor(SC_sampple_id, method = correlation_method)
              # Distances transformation
              SC_sampple_id_dist <- as.dist(SC_sampple_id_corr)
              # Cluster construction of type 'ward.D'
              SC_sampple_id_hc <- hclust(SC_sampple_id_dist, method = hclust_method)
              # Dendogram construction
              SC_sampple_id_dend <- as.dendrogram(SC_sampple_id_hc)
              # Add of colors 
              colors_vector <- colors_vector[order.dendrogram(SC_sampple_id_dend)]
              labels_colors(SC_sampple_id_dend) <- colors_vector
              
              # Clustering plot
              if(ncolors > 15){
                height = ncolors * 0.3
              } else {
                height = 7
              }
              
              # To keep graphic settings
              old_par <- par(no.readonly = TRUE)
              on.exit(suppressWarnings(par(old_par)))
              
              filename <- paste(file_title, ".pdf", sep = "")
              pdf(filename, width = 7, height = height)
              par(cex = 1.2, mar = c(2, 2, 0, 13))
              print(plot(SC_sampple_id_dend, horiz = TRUE))
              dev.off()
              print("Clustering file generated")
            }
            else {
              stop("The third argument must be a character indication the title of the file")
            }
          }
          else{
            stop("Check that the number of colors are equal to the number of samples!")
          }
        }
        else{
          var_options <- paste(colnames(metadata), collapse = ", ")
          stop(paste(c("The second argument must be ONE of the following options:", var_options), collapse = " "))
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