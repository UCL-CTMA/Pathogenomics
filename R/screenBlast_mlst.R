#' Screen within a reference sequence the presence of a querry sequence
#'
#' With this function, the user provides a reference sequence and a querry sequence and the function compute the percentage of the reference which is covered by the querry. This computation is computed using the GenomicRanges Bioconductor objects.
#'
#'
#' @param file_strains The path of the reference sequence that you want to screen. The sequence should be a .fasta file.
#' @param dir_querry the folder containing all fasta files. The sequences are available on : https://pubmlst.org/
#' @param path_profile the path to the file containing the mlst profiles
#' @param path_blastn the path where Blastn is located
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry sequence.
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import dplyr
#' @export screen_Blast_mlst

screen_Blast_mlst <- function (file_strains, dir_querry,path_blastn,path_profile)
{
  profile <- read.table(file = path_profile,header = T)
  querry <- list.files(path = dir_querry,pattern = ".fas",full.names = T)
  querry  <- readDNAStringSet(querry)
  writeXStringSet(x = querry,filepath = "allgenes.fasta")
  querry <- "allgenes.fasta"
  myarg <- paste0(" -subject ",file_strains," -query ",querry," -out blast.txt  -outfmt \"6 qacc qlen length qstart qend pident sacc \"")
  system2(command = path_blastn, args = myarg)
  file.remove("allgenes.fasta")
  blast <- try(read.table("blast.txt"), silent = T)
  file.remove("blast.txt")
  if (class(blast) == "data.frame")
  {
    colnames(blast) <- c( "querry_access", "querry_length", "alignment_lenght", "querry_start", "querry_end","pc_ident", "subject_access")

    blast <- blast %>% filter(pc_ident==100)

    GR <- GRanges(seqnames = blast$querry_access,ranges = IRanges(start =blast$querry_start ,end = blast$querry_end))
    GR_disjoin <- disjoin(GR)

    coverage_df <- tibble(querry_access = seqnames(GR_disjoin)|>as.character(),cover_length = width(GR_disjoin))
    blast = blast %>% select(querry_access,querry_length) %>%
      distinct() %>%
      left_join(coverage_df,by = "querry_access") %>%
      group_by(querry_access,querry_length)%>%
      summarise(cover_length = sum(cover_length))%>%
      mutate(pc_coverage = round(100*cover_length/querry_length,1))%>%
      filter(pc_coverage==100)%>%
      mutate(hit = as.integer(unlist(strsplit(querry_access,"_"))[length(unlist(strsplit(querry_access,"_")))]))%>%
      mutate(querry_access = paste0(unlist(strsplit(querry_access,"_"))[-length(unlist(strsplit(querry_access,"_")))],collapse = "_"))%>%
      select(querry_access,hit)%>%
      pivot_wider(names_from = querry_access, values_from = hit)

    if(length(profile)-1 == length(blast))
    {
      blast <- blast%>%
        left_join(profile,by = names(blast))
    }

    blastvector = as.numeric(blast[1,])
    names(blastvector) = colnames(blast)
    return(blastvector)
  }
  else{warning("error : No hit with blast")
    vec= NA
    names(vec) = "ST"
    return(vec)}
}

