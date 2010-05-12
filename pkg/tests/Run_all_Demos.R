library(RcompHam94)

sourceDir <- function(path, trace = TRUE, ...) {
         for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
            if(trace) cat("#######################\n###",nm,"\n#######################\n")           
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
         }
      }

sourceDir(path="../demo/", trace = TRUE, echo=TRUE) 

unlink("./Rplots.pdf")
