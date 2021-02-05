#' species_fw 
#' 
#' @description Data containing the abundance of 17 species  
#'    expressed in spectral counts. Data generated from an Orbitrap 
#'    Fusion Lumos  Tribrid Mass Spectrometer. The dataset contains 
#'    the metaproteomes from three extraction methods: i) "Q" for Qiagen, 
#'    ii) "FW" for fecal waters, and iii) "Q_FW" for the mixture 
#'    of Qiagen and fecal waters. Data generated in the context of 
#'    the project Microbiome Rapid Access (Université Paris-Saclay).
#' 
#' @docType data
#' 
#' @usage data(species_fw)
#' 
#' @keywords datasets
#' 
#' @format A list of four elements defined as "spectral_count_object" with 
#'    taxonomic annotation generated with the function "crumble_taxonomy"
#'    \describe{
#'        \item{SC_subgroups}{dataframe with 9 samples containing 17 
#'           especies with abundance expressed as spectral counts}
#'        \item{metadata}{information related to the 9 samples from the 
#'           experiment}
#'        \item{peptides_proteins}{information related to each of the 
#'           1557 identified peptides}
#'        \item{type_object}{character indicating the type of object}
#'    }
#' 
"species_fw"