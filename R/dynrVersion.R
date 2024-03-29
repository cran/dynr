
gitVersion <- "v0.1.16-105-g40277ff"

##' Current Version String
##'
##' @param verbose If TRUE, print detailed information to the console (default)
##' 
##' This function returns a string with the current version number of dynr. Optionally (with verbose = TRUE (the default)), it prints a message containing the version of R and the platform.  The primary purpose of the function is for bug reporting.
##' 
##' @return A (length-one) object of class 'package_version'
##' 
##' @examples
##' dynr.version()
##' dynr.version(verbose=FALSE)
##' packageVersion("dynr")
dynr.version <- function (verbose=TRUE) {
	pvers <-  packageVersion("dynr")
	if(verbose){
		msg <- paste("dynr version: ", pvers, " [GIT ", gitVersion, "]", sep="")
		msg <- paste(msg, "\nR version: ", version$version.string, sep="")
		msg <- paste(msg, "\nPlatform: ", version$platform, sep="")
		if ("Darwin" ==Sys.info()["sysname"]){
			msg <- paste(msg, "\nMacOS:", system("sw_vers -productVersion", intern=TRUE))
		}
		message(msg)
	}
	invisible(pvers)
}

