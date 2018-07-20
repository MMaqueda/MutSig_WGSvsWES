
#' Function to select the optimum k value based on different metrics
#' The k value with more votes among the five different metrics is chosen
#' In case there is a tie between different k's, the smallest k is taken
#'
#' @param selCCC    vector of k values selected from CCC
#' @param selRSS    vector of k values selected from RSS
#' @param selFrob   vector of k values selected from Frobenius norm stability
#' @param selSimOriRec    dataframe with the summary of Cos Sim (mean and sd) per k
#' @param nmf_measures    the data frame with the measures from the NMF estimation (several k)
#' @return          A list with two elements: 1st the optimum k and 2nd the votes 
#'                  for each k value
#' 
#' #' Example of usage:
#' kSelection_Consensus(....,nmf_measures=estimateLung$measures)
#'     

kSelection_Consensus <- function(selCCC,selRSS,selFrob,selSimOriRec,nmf_measures){
  
  #Retrieve the best k values based on Cosinus Sim mean and sd
  #Those ones with mean > 0.98
  sel_fromCosSim <- nmf_measures$rank[which(selSimOriRec$mean>0.98)]
  
  #Retrieve the best k values based on Silhoutte.consensus
  #Only those with a value >0.9 are retained 
  sel_fromSilh <- nmf_measures$rank[which(nmf_measures$silhouette.consensus>0.9)]
  
  #Concatenate all values
  all_metrics <- c(selCCC,selRSS,selFrob[[1]],sel_fromCosSim,sel_fromSilh)
  freq_ks <- as.vector(table(all_metrics))
  names(freq_ks) <- paste("Rank",names(table(all_metrics)),sep="")
  
  #Optimum k value: The one with max votes. In case of tie, the minimum k value is taken
  
  return(list(as.integer(names(which.max(table(all_metrics)))),
              freq_ks))
  
}