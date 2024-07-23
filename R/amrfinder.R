#' Title
#'
#' @param organism
#' @param element
#' @param file_strains
#' @param docker_image
#' @param thread
#'
#' @import dplyr
#' @import tidyr
#' @return
#' @export amrfinder
#'
#' @examples
amrfinder <- function(file_strains,organism,element,docker_image = "staphb/ncbi-amrfinderplus:3.12.8-2024-05-02.2",thread = 4){
  myarg <- paste0('run --rm -v ./:/mount_p --cpus=',thread,' ',docker_image,' sh -c "amrfinder -n /mount_p/',file_strains,' -O ', organism,' -d ../amrfinder/data/latest --threads ',thread,' > /mount_p/amrfinder.csv"')
  myarg
  system2(command='docker',args=myarg)

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
