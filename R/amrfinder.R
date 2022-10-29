#' Title
#'
#' @param organism
#' @param path_amrfinder the path where amrfinder is located
#' @param element
#' @param path_database
#' @param file_strains
#'
#' @import dplyr
#' @import tidyr
#' @return
#' @export amrfinder
#'
#' @examples
amrfinder <- function(file_strains,organism,path_amrfinder,element,path_database,thread = 4){
  arg = paste0("-n ",file_strains," -O ", organism," -d ",path_database ," --threads ",thread," > amrfinder.csv")
  system2(command = path_amrfinder,args = arg)
  amr = read.csv(file = "amrfinder.csv",sep = '\t',header = T)
  amr = amr %>%
    select(Gene.symbol,Element.subtype,Class) %>%
    filter(Element.subtype ==element) %>%
    mutate(name = paste0(Class,"_", Gene.symbol))

  amr_vector = rep(1,length(amr$name))
  names(amr_vector) = amr$name
  file.remove("amrfinder.csv")
  return(amr_vector)
}
