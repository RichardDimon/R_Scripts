#' Calculates Nei diversity, given a set of genotypes (gt)
#' and a vector of weight values (equal length to the 
#' number of individuals 
#' 
#' @param gt  - genotype matrix (individuals, loci)
#' @param w   - vector of weights
#' @return Nei diversity, a float value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

common_allele_count <- function(gt, w=NULL, cthresh=1) {

   # estimate allele frequencies
   if (is.null(w)) {

      alt_counts  <- colSums(gt) 
      ref_counts  <- colSums(2-gt)

   } else {

      if (nrow(gt) == length(w)) {
         wm <- matrix(rep(w,ncol(gt)),ncol=ncol(gt),byrow=FALSE) 
         alt_counts  <- colSums(gt*wm, na.rm = TRUE) 
         ref_counts  <- colSums((2-gt)*wm, na.rm = TRUE)

      } else {
         cat("   weights supplied: must have length equal to number of rows in gt \n")
      }
   }

   min_counts <- alt_counts
   for ( i in 1:ncol(gt) ) {

      if ( alt_counts[i] > ref_counts[i] ) {

          min_counts[i] <- ref_counts[i]

      } 
 
   }

   minor_allele_counts <- min_counts
   number_common_alleles <- length(which( min_counts >= cthresh ))

   return( list(number_common_alleles=number_common_alleles, minor_allele_counts=minor_allele_counts) )
}

