#' Title
#'
#'
#' @param file_vfdb
#'
#' @import Biostrings dplyr
#' @return
#' @export
#'
#' @examples
species_list <- function(file_vfdb)
{
  sequences <- readDNAStringSet(file_vfdb)
  nom = names(sequences)
  data = data.frame(nom)
  data <- data %>%
    rowwise()%>%
    mutate(species = unlist(strsplit(nom, '\\['))[length( unlist(strsplit(nom, '\\[')))]) %>%
    mutate(species = paste0(unlist(strsplit(species, ' '))[1:2],collapse = " "))
  list <- gsub(x = unique(data$species),pattern = "\\]",replacement = '')
  return(list)
}

