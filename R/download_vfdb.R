#' Title
#'
#' @param dir_out
#'
#' @import R.utils
#' @return
#' @export
#' @examples
download_vfdb <- function(dir_out,dataset= "core")
{
  if (dataset =="core") {
    download.file(url = "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",
                  destfile = paste0(dir_out,"/VFDB_setA_nt.fas.gz"))
    gunzip(filename = paste0(dir_out,"/VFDB_setA_nt.fas.gz"), remove=T)
  }else if (dataset =="full"){  download.file(url = "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz",
                            destfile = paste0(dir_out,"/VFDB_setA_nt.fas.gz"))
    gunzip(filename = paste0(dir_out,"/VFDB_setA_nt.fas.gz"), remove=T)
  }else{message('You have to choice between "core" or "full" dataset in VFDB database')
  }
}
