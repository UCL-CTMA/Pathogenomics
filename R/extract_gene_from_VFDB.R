#' Title
#'
#' @param file_vfdb
#' @param dir_out
#' @param species
#'
#' @import Biostrings dplyr
#' @return
#' @export
#'
#' @examples
extract_gene <- function(species,file_vfdb,dir_out)
{
  sequences <- readDNAStringSet(file_vfdb)
  sequences <- sequences[grepl(x = names(sequences),pattern = species)]
  nom = names(sequences)
  data <- data.frame(nom)
  data <- data %>%
    rowwise()%>%
    mutate(gene = unlist(strsplit(nom, ' '))[2]) %>%
    mutate(gene = gsub(pattern = '\\(',replacement = "",x = gene))%>%
    mutate(gene = gsub(pattern = '\\)',replacement = "",x = gene))%>%
    mutate(gene = gsub(pattern = '\\/',replacement = "-",x = gene))
  gene = unique(data$gene)
  names(sequences) = data$gene
  allsequences=DNAStringSet()
  for (currentgene in gene)
  {
    currentsequences = sequences[names(sequences)==currentgene]
    if(length(currentsequences)>1){
      for (i in 1:length(currentsequences)) {
        names(currentsequences)[i] = paste0(names(currentsequences)[i],i)
      }

    }
    allsequences = DNAStringSet(c(allsequences,currentsequences))
  }
  path <- paste0(species,"_genes.fasta")
  path <- sub(pattern = " ",replacement = "_",x = path)
  writeXStringSet(x = allsequences,filepath = paste0(dir_out,"/",path))
  message(paste0(dir_out,"/",path))
}

