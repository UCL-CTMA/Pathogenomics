#' Screen within a reference sequence the presence of a querry sequence
#'
#' With this function, the user provides a reference sequence and a querry sequence and the function compute the percentage of the reference which is covered by the querry. This computation is computed using the GenomicRanges Bioconductor objects.
#'
#'
#' @param file_strains The path of the reference sequence that you want to screen. The sequence should be a .fasta file including one or several sequences
#' @param file_querry The path of the querry sequence. he sequence should be a .fasta file including one or several sequences
#' @param docker_image
#' @param pc_id_treshold the minimum percentage of identity obtained with Blast that should be obtained in order to consider that the query match the reference sequence
#' @param thread
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry sequence.
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @import dplyr
#' @export

screen_Blast <- function (file_strains, file_querry,pc_id_treshold,docker_image = "staphb/blast:2.15.0",thread=4)
{
  myarg <- paste0('run --rm -v ./:/mount_p --cpus=',thread,' ',docker_image,' sh -c "blastn -subject /mount_p/',file_strains,' -query /mount_p/',file_querry,' -out /mount_p/blast.txt -outfmt \\"6 qacc qlen length qstart qend pident sacc \\""')
  system2(command='docker',args=myarg)

  blast <- try(read.table("blast.txt"), silent = T)
  file.remove("blast.txt")
  if (class(blast) == "data.frame")
  {
    colnames(blast) <- c( "querry_access", "querry_length", "alignment_lenght", "querry_start", "querry_end","pc_ident", "subject_access")

    blast <- blast %>% filter(pc_ident>pc_id_treshold)

    GR <- GRanges(seqnames = blast$querry_access,ranges = IRanges(start =blast$querry_start ,end = blast$querry_end))
    GR_reduce <- reduce(GR)

    coverage_df <- tibble(querry_access = seqnames(GR_reduce)|>as.character(),cover_length = width(GR_reduce))
    blast = blast %>% select(querry_access,querry_length) %>%
      distinct() %>%
      left_join(coverage_df,by ="querry_access") %>%
      group_by(querry_access,querry_length)%>%
      summarise(cover_length = sum(cover_length))%>%
      mutate(pc_coverage = round(100*cover_length/querry_length,1))


    pc_coverage = blast$pc_coverage
    names(pc_coverage) = blast$querry_access

    return(pc_coverage)
  }
  else{warning("warning : No hit with blast")
    vector = c(1)
    names(vector) = "nohit"
    return(vector)}
}
