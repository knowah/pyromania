
read.pyro.data <- function(file) {
  pyro.data <- tryCatch({
    # try to read with UTF-16 encoding as output by Windows pyromark software
    read.table(file, sep=";", fileEncoding="UTF16", skip=2, header=TRUE, stringsAsFactors=FALSE)
  }, error = function(error_condition) {
    # if that fails, try with default encoding
    try(read.table(file, sep=";", skip=2, header=TRUE, stringsAsFactors=FALSE))
  })
  pyro.data$Well <- as.character(pyro.data$Well)
  
  colnames(pyro.data) <- gsub("\\.$", "", gsub("[\\.]+", "\\.", colnames(pyro.data)))
  
  N.CG <- max(as.numeric(gsub("^Pos\\.([0-9]+)\\.Name", "\\1", 
                              grep("^Pos\\.([0-9]+)\\.Name", colnames(pyro.data), value=TRUE)
  )))
  
  keep.cols <- c(1,2,3, 6 * rep(1:N.CG, each=2) + c(2,3))
  pyro.data <- pyro.data[,keep.cols]
  pyro.data$short.name <- substr(as.character(pyro.data$Assay), 1, 5)
  
  pyro.data
}
