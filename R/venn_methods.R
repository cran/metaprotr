#' venn methods
#' 
#' @description Data containing the subgroups for each logical section of 
#'    the Venn diagram (specific and intersections) from two methods of extraction. 
#'    The extraction methods are: i) "S" for Qiagen, and ii) "S_EF" for
#'    the mixture of Qiagen and fecal waters.
#' 
#' @docType data
#' 
#' @usage data(venn_methods)
#' 
#' @keywords datasets
#' 
#' @format A list of six elements defined as "venn_lists_object", 
#'    generated with the function "plot_venn"
#'    \describe{
#'        \item{S_EF}{ 385 subgroups found in the method 'S_EF'}
#'        \item{S}{ 242 subgroups found the method 'S'}
#'        \item{intersection}{ 201 subgroups in common between the methods 'S' and 'S_EF'}
#'        \item{Specific in S_EF}{ 184 specific subgroups found in the method 'S_EF'}
#'        \item{S}{ 41 specific subgroups found in the method 'S'}
#'        \item{type_object}{character indicating the type of object}
#'    }
#' 
"venn_methods"