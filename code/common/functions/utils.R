suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))


read_header <- function(path, prefix){
  hea <- list()
  if (grepl(".gz$", path)){
    con <- gzfile(path, "rt")
  } else {
    con <- file(path, "r")
  }
  i <- 1
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 | !grepl(paste0("^", prefix), line)) {
      break
    } else {
      hea[[i]] <- line
      i <- i+1
    }
  }
  close(con)

  hea
}


load_table <- function(path, header_prefix=NULL, ...){
  args <- list(...)

  if (!"delim" %in% names(args)){
    if (grepl(".csv", path)){
      args$delim <- ","
    } else if (grepl(".tsv", path) | grepl(".txt", path)){
      args$delim <- "\t"
    }
  }

  if (!is.null(header_prefix)){
    header <- read_header(path, prefix=header_prefix)
    args$skip <- length(header)
  }

  if (grepl(".xlsx$", path)){
    df <- do.call(readxl::read_excel, c(path, args))
  } else {
    if (!"progress" %in% names(args)){
      args$progress <- F
    }
    if (!"show_col_types" %in% names(args)){
      args$show_col_types <- F
    }

    if(grepl(".gz$", path)) {
      file <- base::gzfile(path)
      df <- do.call(readr::read_delim, c(list(file=file), args))
    } else {
      df <- do.call(readr::read_delim, c(path, args))
    }
  }

  df
}
