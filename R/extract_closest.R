#' Title
#'
#' @param lengthconf
#' @param identconf
#' @param offset
#' @param file_strains
#' @param docker_image
#' @param thread
#' @param file_querry
#'
#' @import Biostrings
#' @import GenomicRanges
#' @import dplyr
#' @return
#' @export extract_closest
#'
#' @examples
extract_closest <- function(file_strains,file_querry,lengthconf = 95, identconf =95,offset=0, docker_image = "staphb/blast:2.15.0",thread= 4)
{
  myarg <- paste0('run --rm -v ./:/mount_p --cpus=',thread,' ',docker_image,' sh -c "blastn -subject /mount_p/',file_strains,' -query /mount_p/',file_querry,' -out /mount_p/blast.txt -outfmt \\"7 qacc qlen length pident qstart qend sacc sstart send \\""')
  system2(command='docker',args=myarg)


  blast <- try(read.table('blast.txt', comment.char = '#'),silent=T)
  file.remove("blast.txt")
  if(class(blast)=='data.frame')
  {
    colnames(blast) <- c('querry_access','querry_length','alignment_lenght','pc_ident','querry_start','querry_end','subject_access','subject_start','subject_end')
    blast <- blast %>%
      filter(alignment_lenght>=lengthconf/100*blast$querry_length) %>%
      filter(pc_ident>=identconf) %>%
      arrange(desc(alignment_lenght))

    # GR <- GRanges(seqnames = blast$querry_access,ranges = IRanges(start =blast$querry_start ,end = blast$querry_end),pc_ident = blast$pc_ident)
    # GR_disjoin <- disjoin(GR)
    if(dim(blast)[1]>=1)
    {
      sequence <- readDNAStringSet(file_strains)
      sequence <- sequence[names(sequence)==as.character(blast$subject_access[1])]


      if(blast$subject_start[1]<blast$subject_end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject_start[1]-offset),end=min(blast$subject_end[1]+offset,width(sequence)))}

      if(blast$subject_start[1]>blast$subject_end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject_end[1]-offset),end=min(blast$subject_start[1]+offset,width(sequence)))
      sequence <- reverseComplement(sequence)}
    }
    else{sequence <- ''}
  }
  else{sequence <- ''}
  sequence <- DNAStringSet(sequence)
  return(sequence)
}
