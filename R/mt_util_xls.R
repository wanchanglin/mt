## wll-25-11-2015: This file includes functions for XLS. Most of them are
##   hacked from other packages, especially from 'WriteXLS'. All funstions
##   are internal functionsi,i.e no documents for all these functions.


## ========================================================================
## lwc-17-05-2011: Write list of data frame or matrix to Excel's XLS.
list2xls <- function(x, filename="tmp.xls",FreezeRow=1,row.names = TRUE) {
  tmp <- names(x)
  for (i in tmp)  assign(i, x[[i]])
  mt:::WriteXLS(tmp, ExcelFileName = filename, FreezeRow = FreezeRow,
                row.names = row.names)
  rm(tmp, list=tmp)
  ## Remove tmp and variables whose names are included in tmp
  invisible()
}

## =========================================================================
## lwc-27-05-2010: convert csv to xls.
## lwc-06-04-2011: major change for new version of package WriteXLS.
## Note: This function is hacked from WriteXLS from package WriteXLS.
##       (2.1.0, Date:2010-09-18). For old version, see csv2xls.1
## wll-19-11-2015: '.path.package' is defunct so use 'path.package' instead.
## Arguments:
##  x - A character vector containing the names of one or more CSVs to be
##      converted to the XLS file.
## Note: will remove the original CSV files.
## TO-DO: 1. Check to be sure that each 'x' exist in csv.dir.
##        2. Check csv.dir whether exists or not
## Internal function
## =========================================================================
csv2xls <- function(x, csv.dir, ExcelFileName = "R.xls", perl = "perl",
                    verbose = FALSE, Encoding = c("UTF-8", "latin1"),
                    AdjWidth = FALSE, AutoFilter = FALSE,
                    BoldHeaderRow = FALSE,
                    FreezeRow = 0, FreezeCol = 0) {
  require("WriteXLS", quietly = TRUE)
  Encoding <- match.arg(Encoding)

  SheetNames <- x

  ## Get path to WriteXLS.pl
  Perl.Path <- file.path(path.package("WriteXLS"), "Perl")
  ## Perl.Path <- file.path(.path.package("WriteXLS"), "Perl")
  ## lwc-19-11-2015: '.path.package' is defunct
  Fn.Path <- file.path(Perl.Path, "WriteXLS.pl")

  ## -----------------------------------------------------------------------
  ## Write 'x' (character vector of csv files names) to file
  ## appending csv.dir and ".csv" to each x
  x <- paste(csv.dir, "/", x, ".csv", sep = "")
  write(as.matrix(x), file = paste(csv.dir, "/FileNames.txt", sep = ""))

	## ----------------------------------------------------------------------
  if (verbose) cat("Creating SheetNames.txt\n")
  write(as.matrix(SheetNames),
        file = paste(csv.dir, "/SheetNames.txt", sep = ""))
  if (verbose) cat("\n")

  ## Call Perl script
  cmd <- paste(perl,
               " -I", shQuote(Perl.Path),
               " ", shQuote(Fn.Path),
               " --CSVPath ", shQuote(csv.dir),
               " --verbose ", verbose,
               " --AdjWidth ", AdjWidth,
               " --AutoFilter ", AutoFilter,
               " --BoldHeaderRow ", BoldHeaderRow,
               " --FreezeRow ", FreezeRow,
               " --FreezeCol ", FreezeCol,
               " --Encoding ", Encoding,
               " ", shQuote(ExcelFileName), sep = "")

  ## Call the external Perl script and get the result of the call
  Result <- system(cmd)

  unlink(paste(csv.dir, "/FileNames.txt", sep = ""))
  unlink(paste(csv.dir, "/SheetNames.txt", sep = ""))
  unlink(x)

  if (Result != 0) {
    message("The Perl script 'WriteXLS.pl' failed to run successfully.")
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}

## =========================================================================
## lwc-27-05-2010: convert csv to xls.
## Note: This function is hacked from WriteXLS from package WriteXLS.
##       It is only for package WriteXLS ver 1.9.0(2010-03-20)
## Arguments:
##  x - A character vector containing the names of one or more CSVs to be
##      converted to the XLS file.
## Note: will remove the original CSV files.
## TO-DO: 1. Check to be sure that each 'x' exist in csv.dir.
##        2. Check csv.dir whether exists or not
## Internal function
## ========================================================================
csv2xls.1 <- function(x, csv.dir, ExcelFileName = "R.xls", perl = "perl",
                      verbose = FALSE, Encoding = c("UTF-8", "latin1"),
                      AdjWidth = FALSE, AutoFilter = FALSE,
                      BoldHeaderRow = FALSE,
                      FreezeRow = 0, FreezeCol = 0) {

  require("WriteXLS", quietly = TRUE)
  Encoding <- match.arg(Encoding)

  ## Get path to WriteXLS.pl
  Perl.Path <- file.path(path.package("WriteXLS"), "Perl")
  Fn.Path <- file.path(Perl.Path, "WriteXLS.pl")

  ## -----------------------------------------------------------------------
  ## Write 'x' (character vector of csv files names) to file
  ## appending csv.dir and ".csv" to each x
  x <- paste(csv.dir, "/", x, ".csv", sep = "")
  write(as.matrix(x), file = paste(csv.dir, "/FileNames.txt", sep = ""))

  SN <- FALSE    ## SheetNames. Use default.
  if (verbose) cat("\n")

  ## Call Perl script
  cmd <- paste(perl,
               " -I", shQuote(Perl.Path),
               " ", shQuote(Fn.Path),
               " --CSVPath ", shQuote(csv.dir),
               " --verbose ", verbose,
               " --SN ", SN,
               " --AdjWidth ", AdjWidth,
               " --AutoFilter ", AutoFilter,
               " --BoldHeaderRow ", BoldHeaderRow,
               " --FreezeRow ", FreezeRow,
               " --FreezeCol ", FreezeCol,
               " --Encoding ", Encoding,
               " ", shQuote(ExcelFileName), sep = "")

  ## Call the external Perl script and get the result of the call
  Result <- system(cmd)

  unlink(paste(csv.dir, "/FileNames.txt", sep = ""))
  unlink(x)

  if (Result != 0) {
    message("The Perl script 'WriteXLS.pl' failed to run successfully.")
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}

## ========================================================================
## lwc-06-04-2011: Slight modification of WriteXLS.R (2.1.0,
## Date:2010-09-18) Internal function Write R data frames to an Excel binary
## file using a Perl script
## ========================================================================
WriteXLS <- function(x, ExcelFileName = "R.xls", SheetNames = NULL,
                     perl = "perl", verbose = FALSE,
                     Encoding = c("UTF-8", "latin1"), row.names = FALSE,
                     AdjWidth = FALSE, AutoFilter = FALSE,
                     BoldHeaderRow = FALSE, FreezeRow = 0, FreezeCol = 0,
                     envir = parent.frame()) {

  require("WriteXLS", quietly = TRUE) ## added by LWC

  ## If 'x' is a single name, it is either a single data frame or a list of
  ## data frames. If 'x' is >1 names in a character vector, it is presumed to
  ## be a vector of data frame names. If not a list name, create a list of
  ## data frames from the vector, for consistency in subsequent processing.
  if (length(x) == 1)  {
    TMP <- get(as.character(x), envir = envir)

    ## is TMP a list and not single data frame
    if ((is.list(TMP)) & (!is.data.frame(TMP)))  {
      DF.LIST <- TMP
    } else {
      DF.LIST <- list(TMP)
      names(DF.LIST) <- x
    }
  } else {
    DF.LIST <- sapply(as.character(x), function(x) get(x, envir = envir),
                      simplify = FALSE)
    names(DF.LIST) <- x
  }

  ## Check to be sure that each element of DF.LIST is a data frame

  ## lwc-06-04-2010: comment the following line for allowing matrix format.
  ## if (!all(sapply(DF.LIST, is.data.frame))) stop("One or more of the
  ## objects named in 'x' is not a data frame or does not exist")

  ## Additional checks for Excel 2003 limitations
  ## 256 columns, including rownames, if included
  ## 65,536 rows (including header row)
  if (!all(sapply(DF.LIST, function(x) (nrow(x) <= 65535) & (ncol(x) <= 256))))
    stop("One or more of the data frames named in 'x' exceeds 65535 rows or 256 columns")

  Encoding <- match.arg(Encoding)

  ## -------------------------------------------------------------------------
  ##  Check to see if SheetNames is specified and if so: check for
  ##  duplications they are same length as the number of dataframes CHECK
  ##  TO see if any SheetNames are >31 chars, which is the Excel Limit
  ##  check for invalid characters: []:*?/\ ELSE check to see if first 31
  ##  characters of data frame names are unique
  if (!is.null(SheetNames)) {
    if (any(duplicated(SheetNames))) {
      message("At least one entry in 'SheetNames' is duplicated. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (length(DF.LIST) != length(SheetNames)) {
      message("The number of 'SheetNames' specified does not equal the number of data frames in 'x'")
      return(invisible(FALSE))
    }

    if (any(nchar(SheetNames) > 31)) {
      message("At least one of 'SheetNames' is > 31 characters, which is the Excel limit")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", SheetNames))) {
      message("Invalid characters found in at least one entry in 'SheetNames'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
    names(DF.LIST) <- SheetNames
  } else {
    if (any(duplicated(substr(names(DF.LIST), 1, 31)))) {
      message("At least one data frame name in 'x' is duplicated up to the first 31 characters. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", names(DF.LIST)))) {
      message("Invalid characters found in at least one data frame name in 'x'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
  }

  ## Get path to WriteXLS.pl
  Perl.Path <- file.path(path.package("WriteXLS"), "Perl")
  Fn.Path <- file.path(Perl.Path, "WriteXLS.pl")

  ## Get path for Tmp.Dir for CSV files
  Tmp.Dir <- file.path(tempdir(), "WriteXLS")

  ## Remove Tmp.Dir and Files
  clean.up <- function() {
    if (verbose)
      cat("Cleaning Up Temporary Files and Directory\n\n")
    unlink(Tmp.Dir, recursive = TRUE)
  }

  ## Clean up on function exit
  on.exit(clean.up())

  ## Cleanup now, in case Tmp.Dir still exists from a prior run
  if (file.exists(Tmp.Dir)) {
    if (verbose)
      cat("Cleaning Up Temporary Files and Directory From Prior Run\n\n")

    unlink(Tmp.Dir, recursive = TRUE)
  }

  ## Create Tmp.Dir for new run
  if (verbose)
    cat("Creating Temporary Directory for CSV Files: ", Tmp.Dir, "\n\n")

  dir.create(Tmp.Dir, recursive = TRUE)

  ## ---------------------------------------------------------------------
  ##  Write Comma Delimited CSV files
  for (i in seq(along = DF.LIST)) {
    if (verbose)
      cat("Creating CSV File: ", i, ".csv", "\n", sep = "")

    if (row.names) {
      write.table(DF.LIST[[i]],
                  file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
                  sep = ",", quote = TRUE, na = "",
                  row.names = TRUE, col.names = NA)
    } else {
      write.table(DF.LIST[[i]],
                  file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
                  sep = ",", quote = TRUE, na = "", row.names = FALSE)
    }
  }

  ## Write 'x' (character vector of data frame names) to file
  ## appending Tmp.Dir and ".csv" to each x
  x <- paste(Tmp.Dir, "/", seq(length(DF.LIST)), ".csv", sep = "")
  write(as.matrix(x), file = paste(Tmp.Dir, "/FileNames.txt", sep = ""))

  ## ---------------------------------------------------------------------
  if (verbose)
    cat("Creating SheetNames.txt\n")
  write(as.matrix(names(DF.LIST)),
        file = paste(Tmp.Dir, "/SheetNames.txt", sep = ""))

  if (verbose)
    cat("\n")

  ## Call Perl script
  cmd <- paste(perl,
               " -I", shQuote(Perl.Path),
               " ", shQuote(Fn.Path),
               " --CSVPath ", shQuote(Tmp.Dir),
               " --verbose ", verbose,
               " --AdjWidth ", AdjWidth,
               " --AutoFilter ", AutoFilter,
               " --BoldHeaderRow ", BoldHeaderRow,
               " --FreezeRow ", FreezeRow,
               " --FreezeCol ", FreezeCol,
               " --Encoding ", Encoding,
               " ", shQuote(ExcelFileName), sep = "")

  ## Call the external Perl script and get the result of the call
  Result <- system(cmd)

  ## Check to see if Result != 0 in the case of the failure of the Perl
  ## script This should also raise an error for R CMD check for package
  ## testing on R-Forge and CRAN
  if (Result != 0) {
    message("The Perl script 'WriteXLS.pl' failed to run successfully.")
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}

## ========================================================================
## lwc-21-21-10-2010: Slight modification of WriteXLS.R
## Note: This modification is only for 1.9.0,2010-03-20. The reason for the
##       modification is to allow matraix and save colname and rowname in
##       the Excel file.
##       For the new version, see WriteXLS.R
## Internal function
## ========================================================================
WriteXLS.1 <- function(x, ExcelFileName = "R.xls", SheetNames = NULL,
                       perl = "perl", verbose = FALSE,
                       Encoding = c("UTF-8", "latin1"),
                       AdjWidth = FALSE, AutoFilter = FALSE,
                       BoldHeaderRow = FALSE,
                       FreezeRow = 0, FreezeCol = 0, envir = parent.frame()) {
  ## Check to be sure that each 'x' is a data frame
  ## lwc-21-10-2010: comment the following line.
  ## if (!all(sapply(x, function(i) is.data.frame(get(as.character(i), envir = envir)))))
  ##   stop("One or more of the objects named in 'x' is not a data frame or does not exist")

  require("WriteXLS", quietly = TRUE)
  Encoding <- match.arg(Encoding)

  ## Check to see if SheetNames is specified and if so:
  ##  check for duplications
  ##  they are same length as the number of dataframes
  ##  check to see if any SheetNames are >31 chars, which is the Excel Limit
  ##  check for invalid characters: []:*?/\
  ## ELSE
  ##  check to see if first 31 characters of data frame names are unique
  if (!is.null(SheetNames)) {
    if (any(duplicated(SheetNames))) {
      message("At least one entry in 'SheetNames' is duplicated. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (length(x) != length(SheetNames)){
      message("The number of 'SheetNames' does not equal the number of data frames in 'x'")
      return(invisible(FALSE))
    }

    if (any(nchar(SheetNames) > 31)){
      message("At least one of 'SheetNames' is > 31 characters, which is the Excel limit")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", SheetNames))) {
      message("Invalid characters found in at least one entry in 'SheetNames'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
  } else {
    if (any(duplicated(substr(x, 1, 31)))) {
      message("At least one data frame entry in 'x' is duplicated up to the first 31 characters. Excel worksheets must have unique names.")
      return(invisible(FALSE))
    }

    if (any(grep("\\[|\\]|\\*|\\?|:|/|\\\\", x))) {
      message("Invalid characters found in at least one data frame entry in 'x'. Invalid characters are: []:*?/\\")
      return(invisible(FALSE))
    }
  }

  ## Get path to WriteXLS.pl
  Perl.Path <- file.path(path.package("WriteXLS"), "Perl")
  Fn.Path <- file.path(Perl.Path, "WriteXLS.pl")

  ## Get path for Tmp.Dir for CSV files
  Tmp.Dir <- file.path(tempdir(), "WriteXLS")

  ## Remove Tmp.Dir and Files
  clean.up <- function() {
    if (verbose)
      cat("Cleaning Up Temporary Files and Directory\n\n")

    unlink(Tmp.Dir, recursive = TRUE)
  }

  ## Clean up on function exit
  on.exit(clean.up())

  ## Cleanup now, in case Tmp.Dir still exists from a prior run
  if (file.exists(Tmp.Dir)) {
    if (verbose)
      cat("Cleaning Up Temporary Files and Directory From Prior Run\n\n")

    unlink(Tmp.Dir, recursive = TRUE)
  }

  ## Create Tmp.Dir for new run
  if (verbose)
    cat("Creating Temporary Directory for CSV Files: ", Tmp.Dir, "\n\n")

  dir.create(Tmp.Dir, recursive = TRUE)

  ##  Write Comma Delimited CSV files
  for (i in as.character(x)) {
    if (verbose)
      cat("Creating CSV File: ", i, "\n")

    ## lwc-21-10-2010: comment the following line
    ## write.table(get(i, envir = envir), file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
    ##             sep = ",", quote = TRUE, na = "", row.names = FALSE)

    write.table(get(i, envir = envir), file = paste(Tmp.Dir, "/", i, ".csv", sep = ""),
                sep = ",", quote = TRUE, na = "", row.names = T, col.names=NA)
    ## lwc-21-10-2010: keep the colnames and rownames.
  }

  ## Write 'x' (character vector of data frame names) to file
  ## appending Tmp.Dir and ".csv" to each x
  x <- paste(Tmp.Dir, "/", x, ".csv", sep = "")
  write(as.matrix(x), file = paste(Tmp.Dir, "/FileNames.txt", sep = ""))

  if (!is.null(SheetNames)) {
    if (verbose)
      cat("Creating SheetNames.txt\n")

    write(as.matrix(SheetNames), file = paste(Tmp.Dir, "/SheetNames.txt", sep = ""))
    SN <- TRUE
  } else {
    SN <- FALSE
  }

  if (verbose)
    cat("\n")

  ## Call Perl script
  cmd <- paste(perl,
               " -I", shQuote(Perl.Path),
               " ", shQuote(Fn.Path),
               " --CSVPath ", shQuote(Tmp.Dir),
               " --verbose ", verbose,
               " --SN ", SN,
               " --AdjWidth ", AdjWidth,
               " --AutoFilter ", AutoFilter,
               " --BoldHeaderRow ", BoldHeaderRow,
               " --FreezeRow ", FreezeRow,
               " --FreezeCol ", FreezeCol,
               " --Encoding ", Encoding,
               " ", shQuote(ExcelFileName), sep = "")

  ## Call the external Perl script and get the result of the call
  Result <- system(cmd)

  ## Check to see if Result != 0 in the case of the failure of the Perl
  ## script This should also raise an error for R CMD check for package
  ## testing on R-Forge and CRAN
  if (Result != 0) {
    message("The Perl script 'WriteXLS.pl' failed to run successfully.")
    return(invisible(FALSE))
  } else {
    return(invisible(TRUE))
  }
}

## ==================
## TOC on 25-11-2015
## (1). list2xls
## (2). csv2xls
## (3). csv2xls.1
## (4). WriteXLS
## (5).   clean.up
## (6). WriteXLS.1
## (7).   clean.up
