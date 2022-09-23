#' Title
#'
#' @param prefix
#'
#' @import reticulate
#' @import dplyr
#' @return
#' @export
#' @examples
create_envconda <- function(prefix)
{

  conda_create(envname = paste0("./",prefix))
  conda_install(envname = paste0("./",prefix),packages = "conda")
  conda_install(envname = paste0("./",prefix),packages = "blastn",channel="kantorlab")
  conda_install(envname = paste0("./",prefix),packages = "ncbi-amrfinderplus",forge = T,channel = "bioconda")
  dir.create("amrfinder_database")
  system2(command = paste0(prefix,"/bin/amrfinder_update"),args = " -d amrfinder_database")

  message(paste0("path blastn is: ",prefix,"/bin/blastn"))
  message(paste0("path amrfinder is: ",prefix,"/bin/amrfinder"))
  message("path amrfinder database is : amrfinder_database/latest")
}
