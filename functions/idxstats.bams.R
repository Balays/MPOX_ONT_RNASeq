

#' Alignment IDX statistics from a list of .BAM files.
#'
#' @export

idxstats.bams <- function (bamfiles, bamnames=NULL) {

  require(Rsamtools)
  #require(plyr)

  nbam <- length(bamfiles)

  if (is.null(bamnames) )  {
    pattern  <- '.bam'
    bamnames <- gsub('.*\\/', '', bamfiles)
    bamnames <- gsub(pattern, '', bamnames)
  }

  all.stats <- data.frame(NULL)

  for (i in 1:nbam) {
    stats <- Rsamtools::idxstatsBam(bamfiles[i])
    stats$sample <- bamnames[i]
    all.stats <- plyr::rbind.fill(all.stats, stats)
  }

  return(all.stats)
}
