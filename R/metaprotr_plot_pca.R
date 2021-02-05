#' plot_pca
#'
#' Performs a Principal Components Analysis (PCA) from the spectral counts of the entities 
#' (peptides, subgroups, groups or taxonomic elements) in a 
#' "spectral_count_object" with or without taxonomy. PCA decomposition of
#' high dimensional data allows to observe global effects in two dimensions. For more 
#' details of the used function check dudi.pca from \href{https://CRAN.R-project.org/package=ade4}{ade4}.
#' 
#' @param spectral_count_object List described as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts from 
#'    peptides, subgroups, groups or taxonomic levels. The format of this 
#'    object is similar to that generated from the functions "getsc_specific" and 
#'    "crumble_taxonomy". The PCA projections will be applied to these observations.
#'
#' @param colors_var Character indicating the name of one column 
#'    from metadata. The samples will be represented in different colors
#'    in function of the levels of this variable (ex. conditions).
#'
#' @param pc_components Two numeric values indicating two principal 
#'    components to be analyzed.
#' 
#' @param force Logic value set as FALSE by default in order to ask 
#'    permission to create a file in the workstation of the user.
#'
#' @return A pdf file containing the results of PCA applied to the two provided 
#'    principal components. Including a bar plot indicating the percentage of variance
#'    per principal component.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(fecal_waters)
#' plot_pca(fecal_waters, "Methods", c(1, 2))
#' 
#' data(species_fw)
#' plot_pca(species_fw, "Methods", c(1, 3))
#' 
#' data(species_annot_fw)
#' plot_pca(species_annot_fw, "Condition", c(1, 2))
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

plot_pca <-
function(spectral_count_object, colors_var, pc_components, force = FALSE){
  
  variables <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <-
      select.list(
        c("yes", "no"), 
        title = "Do you allow to create a pdf file with PCA results?"
      )
    
    if (response == "yes") {
      if (length(spectral_count_object) == 4) {
        if (spectral_count_object[[4]] == "spectral_count_object") {
          
          colors_var <- colors_var
          metadata <- spectral_count_object[[2]]
          
          if (colors_var %in% colnames(metadata) == TRUE) {
            
            if ((length(pc_components) == 2) & 
                (length(unique(pc_components)) == 2) &
                (all(is.numeric(pc_components)) == TRUE)) {
              
              if((round(pc_components[1], 0) > 0) & 
                 (round(pc_components[2], 0) > 0) &
                 (round(pc_components[2], 0) < length(metadata$SampleID))) {
                
                # SC data 
                sc_data <- spectral_count_object[[1]]
                # colnames(sc_data) <- metadata[["SampleID"]]
                # sc_data <- sapply(split.default(sc_data, names(sc_data)), rowMeans)
                # sc_data <- as.data.frame(sc_data)
                # for aesthetics in plot
                pepsprots <- spectral_count_object[[3]]
                if("species" %in% colnames(pepsprots)){
                  sc_origin <- "Species"
                }
                else if("genus" %in% colnames(pepsprots)){
                  sc_origin <- "Genus"
                }
                else if("family" %in% colnames(pepsprots)){
                  sc_origin <- "Family"
                }
                else if("order" %in% colnames(pepsprots)){
                  sc_origin <- "Order"
                }
                else if("class" %in% colnames(pepsprots)){
                  sc_origin <- "Class"
                }
                else if("phylum" %in% colnames(pepsprots)){
                  sc_origin <- "Phylum"
                }
                else if("superkingdom" %in% colnames(pepsprots)){
                  sc_origin <- "Superkingdom"
                }
                else if(names(spectral_count_object[1]) == "SC_groups"){
                  sc_origin <- "Groups"
                }
                else if(names(spectral_count_object[1]) == "SC_subgroups"){
                  sc_origin <- "Subgroups"
                }
                else{
                  sc_origin <- "Peptides"
                }
                
                # PCA
                result_pca <- dudi.pca(t(sc_data), center = TRUE, scale = TRUE, scannf = FALSE, 
                                       nf = length(metadata$SC_name) - 1
                )
                variance <- round((result_pca$eig / sum(result_pca$eig) * 100), 2)
                pca_vars <- result_pca$l1
                pca_vars$SC_name <- row.names(pca_vars)
                pca_vars <- merge(pca_vars, 
                                  metadata, 
                                  by = "SC_name", 
                                  all.x = TRUE
                )
                pca_vars[["SC_name"]] <- as.factor(pca_vars[["SC_name"]])
                pca_vars$variables <- pca_vars[[colors_var]]
                xdata <- paste(c("RS", pc_components[1]), collapse = "")
                ydata <- paste(c("RS", pc_components[2]), collapse = "")
                pca_vars$xdata <- pca_vars[[xdata]]
                pca_vars$ydata <- pca_vars[[ydata]]
                variance1 <- variance[pc_components[1]]
                variance2 <- variance[pc_components[2]]
                title <- paste(
                  c("PCA from the sum of spectral counts per", sc_origin),
                  collapse = " "
                )
                
                # To keep graphics settings
                old_par <- par(no.readonly = TRUE)
                on.exit(suppressWarnings(par(old_par)))
                
                plot <- ggplot(pca_vars) + 
                  geom_point(aes(
                    x = xdata, y = ydata, color = variables), 
                    alpha = 0.7, size = 4) +
                  labs(x = paste("PC", pc_components[1] , "(", variance1, " %)\n"), 
                       y = paste("PC", pc_components[2], "(", variance2, " %)\n"), 
                       title = title) +
                  # xlim(-1.1, 1.1) +
                  # ylim(-1.1, 1.1) + 
                  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
                  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
                  theme_minimal()
                filename <- paste("PCA", sc_origin, "pc", 
                                  pc_components[1], pc_components[2], sep = "_")
                filename <- paste(filename, ".pdf", sep = "")
                
                # plot size in function of the number of principal components
                if(length(variance) > 15){
                  width = length(variance) * 0.35
                } 
                else {
                  width = 11
                }
                height = width * 0.8
                pdf(filename, width = width, height = height)
                print(plot)
                max_pc <- length(variance) - 1
                bar_coor <- barplot(
                  variance[1:max_pc], ylab = "Variance (%)", 
                  xlab = "Principal components",
                  main = "Contribution of the Principal Components",
                  ylim = c(0, max(variance) + 5),
                  names.arg = 1:max_pc
                )
                text(x = bar_coor, 
                     y = variance[1:max_pc], label = variance[1:max_pc], 
                     cex = 1, pos = 3, col = "blue",
                )
                box()
                dev.off()
                print(paste("Graph", filename, "was created!", sep = " "))
              }
              else{
                max_val <- length(metadata$SampleID) - 1
                stop(
                  paste(
                    c("Numeric values cannot be lower than 0 or higher than", 
                      max_val), collapse = " ")
                )
              }
            }
            else{
              stop("Precise TWO different numeric values on third argument")
            }
          }
          else{
            options <- colnames(metadata)
            stop(paste(c("Change colors_var for ONE of these options: ", 
                         options), collapse = "' '")
            )
          }
        }
        
        else {
          stop("Invalid spectral count object")
        }
      }
      
      else{
        stop("Invalid spectral count object")
      }
    }
    else{
      stop("No file was created")
    }
  }
}
