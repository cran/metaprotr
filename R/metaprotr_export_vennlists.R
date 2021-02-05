#' export_vennlists
#'
#' Exports as csv files the elements (groups, subgroups, peptides or 
#' taxonomic levels) generated from the function "plot_venn".
#' 
#' @param venn_lists_object List defined as "venn_lists_object" containing the
#'    elements (peptides, subgroups, groups or taxonomic elements) 
#'    generated with the function "plot_venn".
#'
#' @param output_repo Character indicating the path of a previously 
#'    created directory where the lists will be exported. 
#'    This parameter is optional.
#'
#' @param force Logic value set at FALSE by default in order 
#'    to ask permission to create csv files in the workstation of the user.
#'
#' @return csv files containing the elements present on each logic section 
#'    (specific and intersections) from the list defined as "venn_lists_object".
#' 
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' 
#' data(venn_methods)
#' 
#' export_vennlists(venn_methods)
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

export_vennlists <- 
  
function(venn_lists_object, output_repo = NULL, force = FALSE){
  
  venn_lists_object <- venn_lists_object
  
  if (!force && interactive()) {
    response <-
      select.list(c("yes", "no"), title = "Do you allow to create csv files containing the Venn diagram information?")
    
    if (response == "yes") {
      
      if (rev(venn_lists_object)[[1]] == "venn_lists_object" &
          rev(names(venn_lists_object))[[1]] == "type_object") {
        
        print("=============================================")
        for (list in 1:(length(venn_lists_object) - 1)) {
          name_cond <- names(venn_lists_object[list])
          list <- venn_lists_object[[list]]
          filename <- paste("Venn_list", name_cond, sep = "_")
          # If repository is provided, csv files are exported there
          if (length(output_repo) == 0) {
            filename <- paste(filename, ".csv", sep = "")
          }
          else{
            filename <- paste(output_repo, "/", filename, ".csv", sep = "")
          }
          write.table(list, filename, append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
          print(paste(filename, "generated", sep = " "))
        }
        print("All files were exported")
        print("=============================================")
      }
      
      else{
        stop("Invalid venn lists object. The input should come from plot_venn function.")
      }
    }
    else{
      stop("No file was created")
    }
  }
}