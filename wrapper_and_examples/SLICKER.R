SLICKER <- function(target, rtime, p1, p2 = NULL, p3 = NULL, p4 = NULL, 
              p5 = NULL, p6 = NULL, p7 = NULL, p8 = NULL, p9 = NULL, 
              p10 = NULL, ...) {
#
#RUN_SLICKER - runs Fortran SLICKER reconstruction code "reconstruct"
#
#   Brief description : This is a R wrapper that generates input files to 
#   run the Fortran SLICKER reconstruction code described in 
#   Roberts et al. (2023).
#
#   User inputs:
#       target     2d array; 1st column is time, second column is the corresponding data value
#       rtime      1d vector of list of times that you want the reconstruction for
#       p1         2d array; 1st column is time, 2nd column is proxy value
#       p[2-10]    [optional] up to an additional 9 input proxy time series
#       ...        [optional] parameters, only needed to change defaults
#
#   Input arrays may be different lengths. Data for each MUST be in strictly increasing time order.
#
#   Optional parameters and default values
#       num_ensemble    Number of ensemble members (default 4096)
#       slick_width     Slick width parameters (default [0.4 1.6])
#       time_limit      Time limit (seconds) for solution of individual 
#                        ensemble member (default 60 seconds)
#       tol             Solution tolerance (default 1e-6)
#       nonlinear       Attempt a non-linear solution (default "n")
#       frac_subset     Fraction of ensemble with best stationarity used in
#                        reported solution (default 0.5)
#       proxy_inversion Test for significantly stronger correlation if 
#                        individual proxies are inverted (default "n" for 
#                        linear solution and "y" for non-linear)
#
#   Output:
#       a m*4 array (where m is the length of "rtime")
#           column 1 is the reconstruction times
#           column 2 is the central tendency of the SLICKER ensemble solution
#           column 3 is the 95% confidence interval of the central tendency
#           column 4 is the spread of the ensemble
#
#   Usage:
#       recon <- SLICKER(target,rtime,proxy1,proxy2)
#
#
#   Examples:
#       recon <- SLICKER(target,rtime,proxy1,num_ensemble=512)
#       recon <- SLICKER(target,rtime,proxy1,proxy2,proxy3,time_limit=20,tol=1e-7)
#
#
#   Citations:
#   If you use this R wrapper, please cite
#   Roberts et al. (submitted), 2023.
#
#   Author Info
#   This function was written by Jason Roberts and John French     
#   May 2023
#

  defaults <- list(...)

  if (!is.null(defaults$num_ensemble)) {
    num_ensemble <- defaults$num_ensemble
  } else {
    num_ensemble <- 4096
  }
  if (!is.null(defaults$slick_width)) {
    slick_width <- defaults$slick_width
  } else {
    slick_width <- c(0.4,1.6)
  }
  if (!is.null(defaults$time_limit)) {
    time_limit <- defaults$time_limit
  } else {
    time_limit <- 60
  }
  if (!is.null(defaults$tol)) {
    tol <- defaults$tol
  } else {
    tol <- 1e-6
  }
  if (!is.null(defaults$nonlinear)) {
    isnonlinear <- tolower(defaults$nonlinear)
  } else {
    isnonlinear <- "n"
  }
  if (!is.null(defaults$frac_subset)) {
    frac_subset <- defaults$frac_subset
  } else {
    frac_subset <- 0.5
  }
  if (!is.null(defaults$proxy_inversion)) {
    ispinversion <- tolower(defaults$proxy_inversion)
  } else if (isnonlinear == "y") {
    ispinversion <- "y"
  } else {
    ispinversion <- "n"
  }

  num_proxy_files <- 1
  proxy <- list()
  proxy[[1]] <- p1
  proxy_list <- list(p2, p3, p4, p5, p6, p7, p8, p9, p10)
  for (i in proxy_list) {
    if (!is.null(i)) {
      num_proxy_files <- num_proxy_files + 1
      proxy[[num_proxy_files]] <- i
    }
  }

  if (length(target) != 2) {
    stop("SLICKER error message: reconstruction target should be Mx2.")
  }

  for (i in 1:num_proxy_files) {
    if (length(proxy[[i]]) != 2) {
      stop("SLICKER error message: reconstruction proxy ", i, " should be Mx2.")
    }
  }

# Filename of reconstruction target data
basename <- paste0("SLICKER_", format(Sys.time(), "%b_%d_%Y_%H_%M_%S"), "_")
rfilename <- paste0(basename, "target.txt")

# Filenames of proxy files
cat_pfilenames <- toString(num_proxy_files)
pfilenames <- list()
for (i in 1:num_proxy_files) {
    pfilenames[[i]] <- paste0(basename, "proxy", toString(i), ".txt")
    cat_pfilenames <- paste(cat_pfilenames, pfilenames[[i]], sep = "\n")
}

# Filename of times to be used for reconstruction
tfilename <- paste0(basename, "timebase.txt")

# Filename of output file
ofilename <- paste0(basename, "output.txt")

# Filename of master input file
mfilename <- paste0(basename, "input_keyboard_file.txt")

# Generate user input file
fid <- file(mfilename)
writeLines(c(rfilename, cat_pfilenames, tfilename, toString(num_ensemble),
           toString(slick_width), toString(time_limit), toString(tol),
           ofilename, isnonlinear, toString(frac_subset), ispinversion), fid)
close(fid)

# Write target file
write.table(target, rfilename, sep = " ", row.names = FALSE, col.names = FALSE)

# Write reconstruction time file
write.table(rtime, tfilename, sep = " ", row.names = FALSE, col.names = FALSE)

# Write proxy files
for (i in 1:num_proxy_files) {
  write.table(proxy[[i]], pfilenames[[i]], sep = " ", row.names = FALSE, 
              col.names = FALSE)
}

# Run Fortran code
sys <- Sys.info()[1]

if (sys=="Windows") {
  print("OS:Windows")
  system_cmd <- paste0("..\\reconstruct_Win.exe < ", mfilename, " > console.txt")
  shell(system_cmd, wait=TRUE)
}
else if (sys=="UNIX" || sys=="Linux") {
  print("OS:UNIX")
  exec_dir <- getSrcDirectory(SLICKER)
  system_cmd <- paste0(exec_dir, "/reconstruct < ", mfilename, " > /dev/null")
  system(system_cmd)
}

# Read the output file
data <- read.delim(ofilename, header = FALSE, sep = "", comment.char = "!",
                   strip.white = TRUE)

return(data)
}


