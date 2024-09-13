get_minor_allele_frequencies <- function( gt ) {

   alt_freq  <- colSums(gt, na.rm = TRUE) / (2*nrow(gt)) #dont use na.rm unless you are sure. number of SNPs will be cut off
   ref_freq  <- colSums(2-gt, na.rm = TRUE) / (2*nrow(gt))

   min_freq <- alt_freq
   for ( i in 1:ncol(gt) ) {

      if ( alt_freq[i] > ref_freq[i] ) {

          min_freq[i] <- ref_freq[i]

      } 
 
   }
   return(min_freq)
} 
# 
# apply(dd, 2, function(x){ length(which(!is.na(x))) })
# 
# dd <- data.frame(x=c(1,2,3,4),y=c(3,4,5,2))
