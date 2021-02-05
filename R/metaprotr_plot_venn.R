#' plot_venn
#'
#' Generates a Venn diagram comparing up to 3 conditions. The lists of 
#' elements for each condition are also returned as a "venn_lists_object".
#' 
#' @param spectral_count_object List defined as "spectral_count_object" 
#'    containing dataframes with abundance expressed as spectral counts 
#'    from peptides, subgroups, groups or taxonomic levels. 
#'    The format of this object is similar to that generated from the 
#'    functions "getsc_specific" and "crumble_taxonomy".
#'
#' @param target_variable Character indicating the name of the explanatory 
#'    variable that contains the conditions to be compared. This value corresponds 
#'    to the name of one column from the metadata dataframe.
#'
#' @param list_conditions Atomic vector indicating the conditions 
#'    to be compared. The provided elements (2 or 3) must be present in the 
#'    variable indicated as "target_variable".
#'
#' @param force Logic value set at FALSE by default in order to ask 
#'    permission to create a pdf file in the workstation of the user. 
#'
#' @return A Venn diagram (pdf) and a list defined as "venn_list_object"
#'    containing the elements (peptides, soubgroups, groups or taxonomic levels) 
#'    for each logical section of the Venn diagram (specific and intersections).
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' data(fecal_waters)
#' venn_QFW1_Q1 <- plot_venn(fecal_waters, "SC_name", c("Q1_FW1", "Q1"))
#' 
#' data(species_fw)
#' venn_all <- plot_venn(species_fw, "Methods", c("S_EF", "S", "EF"))
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

plot_venn <-

function(spectral_count_object, target_variable, list_conditions, force = FALSE){
  
  conditions <- uniques <- tx <- ty <- x <- y <- NULL
  spectral_count_object <- spectral_count_object
  
  if (!force && interactive()) {
    response <-
      select.list(c("yes", "no"), title = "Do you allow to create a pdf file with a Venn diagram?")
    
    if (response == "yes") {
      
      if (length(spectral_count_object) == 4) {
        
        if (spectral_count_object[[4]] == "spectral_count_object") {
          
          target_variable <- target_variable
          metadata <- spectral_count_object[[2]]
          
          if (target_variable %in% colnames(metadata) == TRUE) {
            
            if (all((unique(list_conditions)) %in% metadata[[target_variable]]) == TRUE) {
              
              if ((length(list_conditions) == 2) &
                  (length(unique(list_conditions)) == 2)) {
                
                sc_data <- spectral_count_object[[1]]
                list_conditions <- list_conditions
                
                # change sample names in function of the selected variable
                #colnames(sc_data) <- metadata$SampleID
                colnames(sc_data) <- metadata[[target_variable]]
                
                # aggregate columns with the same name in target_variable
                sc_data <- sapply(split.default(sc_data, names(sc_data)), rowSums)
                sc_data <- as.data.frame(sc_data)
                cond1 <- list_conditions[1]
                cond2 <- list_conditions[2]
                
                # Venn lists
                set1 <- row.names(sc_data[which(sc_data[[cond1]] != 0), ])
                set2 <- row.names(sc_data[which(sc_data[[cond2]] != 0), ])
                
                # Interection
                set_12 <- intersect(set1, set2)
                nb_s12 <- as.character(length(set_12))
                
                # Differences
                set_unique_1 <- setdiff(set1, set2)
                nb_s1 <- as.character(length(set_unique_1))
                
                set_unique_2 <- setdiff(set2, set1)
                nb_s2 <- as.character(length(set_unique_2))
                
                # Plot logic
                # DataFrame to specify circle and text positions
                data <- data.frame(
                  x = c(-1, 0),
                  y = c(1, 1),
                  tx = c(-2, 0.8),
                  ty = c(1, 1),
                  conditions = list_conditions,
                  uniques = c(nb_s1, nb_s2)
                )
                
                # To keep graphics settings
                old_par <- par(no.readonly = TRUE)
                on.exit(suppressWarnings(par(old_par)))
                
                venn_diagram <- ggplot(
                  data, aes(x0 = x, y0 = y, r = 1.3, fill = conditions)) +
                  geom_circle(alpha = 0.4, size = 1, color = "black", show.legend = TRUE) +
                  geom_text(aes(x = tx, y = ty, label = uniques), size = 7) +
                  annotate(geom = "text", x = -0.5, y = 1, label = nb_s12, color = "black", size = 7) +
                  theme_void()
                venn_diagram <- venn_diagram + theme(legend.position = "bottom")
                
                filename <- paste("VennPlot", cond1, "VS", cond2, sep = "_")
                filename <- paste(filename, ".pdf", sep = "")
                pdf(filename, width = 6.5, height = 4.5)
                print(venn_diagram)
                dev.off()
                
                # Return an object with the lists of the Venn diagram
                # OUTPUT DATA 
                output <-
                  list(set1,
                       set2,
                       set_12,
                       set_unique_1,
                       set_unique_2,
                       "venn_lists_object")
                names(output)[1] <- cond1
                names(output)[2] <- cond2
                names(output)[3] <- "intersection"
                names(output)[4] <- paste("Specific in ", cond1, sep = "")
                names(output)[5] <- paste("Specific in ", cond2, sep = "")
                names(output)[6] <- "type_object"
                # Summary message
                print("=====================================================")
                print(paste("Condition", cond1, "has", as.character(length(set1)), "elements", sep = " "))
                print(paste("Condition", cond2, "has", as.character(length(set2)), "elements", sep = " "))
                print(paste("Venn object and", filename, "were created!", sep = " "))
                print("=====================================================")
                on.exit(dev.off(), add = TRUE)
                return(output)
              }
              
              else if ((length(list_conditions) == 3) &
                       (length(unique(list_conditions)) == 3)) {
                sc_data <- spectral_count_object[[1]]
                list_conditions <- list_conditions
                
                # change sample names in function of the selected variable
                #colnames(sc_data) <- metadata$SampleID
                colnames(sc_data) <- metadata[[target_variable]]
                
                # aggregate columns with the same name in target_variable
                sc_data <- sapply(split.default(sc_data, names(sc_data)), rowSums)
                sc_data <- as.data.frame(sc_data)
                cond1 <- list_conditions[1]
                cond2 <- list_conditions[2]
                cond3 <- list_conditions[3]
                
                # Venn lists
                set1 <- row.names(sc_data[which(sc_data[[cond1]] != 0), ])
                set2 <- row.names(sc_data[which(sc_data[[cond2]] != 0), ])
                set3 <- row.names(sc_data[which(sc_data[[cond3]] != 0), ])
                
                # Interections
                set_12 <- intersect(set1, set2)
                set_123 <- intersect(set_12, set3)
                set_12 <- setdiff(set_12, set_123)
                set_13 <- intersect(set1, set3)
                set_13 <- setdiff(set_13, set_123)
                set_23 <- intersect(set2, set3)
                set_23 <- setdiff(set_23, set_123)
                
                nb_s123 <- as.character(length(set_123))
                nb_s12 <- as.character(length(set_12))
                nb_s13 <- as.character(length(set_13))
                nb_s23 <- as.character(length(set_23))
                
                # Differences
                set_unique_1 <- setdiff(set1, set2)
                set_unique_1 <- setdiff(set_unique_1, set3)
                nb_s1 <- as.character(length(set_unique_1))
                
                set_unique_2 <- setdiff(set2, set1)
                set_unique_2 <- setdiff(set_unique_2, set3)
                nb_s2 <- as.character(length(set_unique_2))
                
                set_unique_3 <- setdiff(set3, set1)
                set_unique_3 <- setdiff(set_unique_3, set2)
                nb_s3 <- as.character(length(set_unique_3))
                
                # Plot logic
                # DataFrame to specify circle and text positions
                data <- data.frame(
                  x = c(0, 1,-1),
                  y = c(-0.5, 1, 1),
                  tx = c(0, 1.5,-1.5),
                  ty = c(-1, 1.3, 1.3),
                  conditions = list_conditions,
                  uniques = c(nb_s1, nb_s2, nb_s3)
                )
                
                # To keep graphics settings
                old_par <- par(no.readonly = TRUE)
                on.exit(suppressWarnings(par(old_par)))
                
                venn_diagram <- ggplot(
                  data, aes(x0 = x, y0 = y, r = 1.7, fill = conditions)) +
                  geom_circle(alpha = 0.4, size = 1, color = "black", show.legend = TRUE) +
                  geom_text(aes(x = tx, y = ty, label = uniques), size = 7) +
                  annotate(geom = "text", x = 0, y = 1.5, label = nb_s23, color = "black", size = 5) +
                  annotate(geom = "text", x = -0.9, y = 0, label = nb_s13, color = "black", size = 5) +
                  annotate(geom = "text", x = 0.9, y = 0, label = nb_s12, color = "black", size = 5) +
                  annotate(geom = "text", x = 0, y = 0.5, label = nb_s123, color = "blue", size = 6) +
                  theme_void()
                venn_diagram <- venn_diagram + theme(legend.position = "bottom")
                
                filename <- paste("VennPlot", cond1, "VS", cond2, "VS", cond3, sep = "_")
                filename <- paste(filename, ".pdf", sep = "")
                pdf(filename, width = 6.5, height = 5.5)
                print(venn_diagram)
                dev.off()
                
                # Return an object with the lists of the Venn diagram
                # OUTPUT DATA
                output <- list(
                  set1,
                  set2,
                  set3,
                  set_12,
                  set_13,
                  set_23,
                  set_123,
                  set_unique_1,
                  set_unique_2,
                  set_unique_3,
                  "venn_lists_object"
                )
                names(output)[1] <- cond1
                names(output)[2] <- cond2
                names(output)[3] <- cond3
                names(output)[4] <- paste("Intersection", cond1, cond2, sep = "_")
                names(output)[5] <- paste("Intersection", cond1, cond3, sep = "_")
                names(output)[6] <- paste("Intersection", cond2, cond3, sep = "_")
                names(output)[7] <- paste("Intersection", cond1, cond2, cond3, sep = "_")
                names(output)[8] <- paste("Specific in ", cond1, sep = "")
                names(output)[9] <- paste("Specific in ", cond2, sep = "")
                names(output)[10] <- paste("Specific in ", cond3, sep = "")
                names(output)[11] <- "type_object"
                
                # Summary message
                print("=====================================================")
                print(paste("Condition", cond1, "has", as.character(length(set1)), "elements", sep = " "))
                print(paste("Condition", cond2, "has", as.character(length(set2)), "elements", sep = " "))
                print(paste("Condition", cond3, "has", as.character(length(set3)), "elements", sep = " "))
                print(paste("Venn object and", filename, "were created!", sep = " "))
                print("=====================================================")
                return(output)
              }
              
              else{
                stop("Precise 2 or 3 conditions in list_conditions argument")
              }
            }
            
            else{
              vars <- paste(list_conditions, collapse = ", ")
              stop(paste(c("The variables [", vars, "] must be present in: '", target_variable, "'"), collapse = " "))
            }
          }
          
          else{
            options <- colnames(metadata)
            stop(paste(c("Change target_variable for ONE of these options: ", options), collapse = "' '"))
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
