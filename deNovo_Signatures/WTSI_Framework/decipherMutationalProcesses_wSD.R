decipherMutationalProcesses_wSD <- function (input, params) 
{
  options(warn = -1)
  requiredElements <- c("cancerType", "mutCounts", "sampleNames", 
                        "mutTypes")
  if (!("input" %in% ls()) | !(is.list(input)) | sum(requiredElements %in% 
                                                     names(input)) != length(requiredElements)) 
    stop("Malformed Input / Input does not include all the required fields!")
  input <- input
  if (is.numeric(params$num.processes.toextract)) {
    params$analyticApproach <- "denovo"
  }
  else {
    stop("An error occurred!")
  }
  if (parallel::detectCores() < 4) 
    stop("Weak computational environment. At least 4 cores are required for running this framework...")
  if (params$analyticApproach == "denovo") {
    deconvData <- deconvoluteMutCounts(input.mutCounts = input$mutCounts, 
                                       params = params)
  }
  else {
    stop("An error occurred!")
  }
  mutProcesses <- list()
  mutProcesses$input <- input
  mutProcesses$params <- params
  mutProcesses$allProcesses <- deconvData$Wall
  mutProcesses$allExposures <- deconvData$Hall
  mutProcesses$idx <- deconvData$idx
  mutProcesses$processes <- deconvData$processes
  mutProcesses$processesStd <- deconvData$processesStd
  mutProcesses$exposures <- deconvData$exposure
  mutProcesses$exposuresStd <- deconvData$exposureStd
  mutProcesses$mutCountErrors <- deconvData$mutCountErrors
  mutProcesses$mutCountReconstructed <- deconvData$mutCountReconstructed
  mutProcesses$processStab <- deconvData$processStab
  mutProcesses$processStabAvg <- deconvData$processStabAvg
  options(warn = 0)
  return(mutProcesses)
}