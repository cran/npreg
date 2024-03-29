# Some code (for startup message) is re-purposed from the R package
# "mclust" (Scrucca et al., 2016) https://cran.r-project.org/package=mclust

npregStartupMessage <- 
  function(){
    msg <- c(paste0("                             
  _ __  _ __  _ __ ___  __ _ 
 | '_ \\| '_ \\| '__/ _ \\/ _` |
 | | | | |_) | | |  __/ (_| |
 |_| |_| .__/|_|  \\___|\\__, |
       |_|             |___/ version ", 
                    packageVersion("npreg"), "\n"),
             "\nType 'citation(\"npreg\")' to cite this package.\n")
    return(msg)
  }

.onAttach <- 
  function(lib, pkg){
    msg <- npregStartupMessage()
    if(!interactive()) msg[1] <- paste("Package 'npreg' version", packageVersion("npreg"))
    packageStartupMessage(msg)      
    invisible()
  }




