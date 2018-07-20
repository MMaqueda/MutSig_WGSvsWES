#Function copy-paste from SomaticSignatures package
#Encountered problems with the package installation so...direct solution

normalizeMotifs <- function (x, norms) 
{
  s = BStringSet(rownames(x))
  base_motif = subseq(s, 4, 6)
  subseq(base_motif, 2, 2) = subseq(s, 1, 1)
  bs = as(base_motif, "character")
  stopifnot(all(bs %in% names(norms)))
  idx = match(bs, names(norms))
  sss = as.vector(norms[idx])
  stopifnot(nrow(x) == length(sss))
  y = x * sss
  return(y)
}