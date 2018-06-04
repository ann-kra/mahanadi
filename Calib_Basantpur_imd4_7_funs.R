# mops_ppse functions Basantpur 2017-09-19

Sys.setenv(TZ="UTC")
library("mops")
library("ppso")

# goodness of fit ----
# modified nash-sutcliffe efficiency

gof_func <- function(obs, sim) {
  NSE = 1 - mean((sim-obs)^2)/var(obs)
  #NSE2 = (-1. + mean((obs - sim)^2) / var(obs) )
  pBias = (sum(sim-obs)/sum(obs))*100
  #mNSE = (NSE * ((100-abs(pBias))/100))*(-1)    # modified NS
  rmse= sqrt(mean((sim-obs)^2))
  ifelse(NSE>0, mNSE <- (NSE * ((100-abs(pBias))/100))*(-1),
                mNSE <- (NSE * ((100+abs(pBias))/100))*(-1))
  all = c(mNSE=mNSE, pBias=pBias, rmse=rmse, NSE=NSE)
  return(all)
}



# computes total error from Q-error and ETR-error
total_error <- function(error1, error2) {
  # just reduce mNSE for Q by pBias
  mNSE_Q = error1["mNSE"]
  pBias_ETR = error2["pBias"]
  # ATTENTION: mNSE is already negative here - the more negative the better!
  if (mNSE_Q < 0) {
    total = mNSE_Q * (100. - abs(pBias_ETR)) / 100.
  } else {
    total = mNSE_Q * (100. + abs(pBias_ETR)) / 100.
  }
  return(total)
}



# Miscellaneous stuff ----

# Process arguments
updatePars_byTemplate <- function(parameters, my_model_args) {
  stopifnot(all(c("file_template","file_result","char_open","char_close")
                %in% names(my_model_args)))
  n= mops::update_template(file_template= my_model_args["file_template"], 
                           vect_updates = parameters,
                           char_open    = my_model_args["char_open"], 
                           char_close   = my_model_args["char_close"],
                           file_result  = my_model_args["file_result"], 
                           overwrite    = TRUE)
  if(n==0) stop(paste0("Updating of '",my_model_args["file_template"],"' failed."))
  return(0)
}


# remove temporary outputs
empty_outputfolder <- function(...) {
  
  if (!dir.exists(outdir)) stop("Directory '", outdir, "' does not exist.")
  files <- list.files(path=outdir, full.names=TRUE)
  tlh <- tools::file_ext(files) %in% c("txt", "html", "log")
  files <- files[tlh] # keep all other files, for safety
  file.remove(files)
  return(0)  # debug
}


# dummy first/final functions to be used whenever needed
dummygc <- function(...) {gc(verbose=FALSE) ;  return(0)}


# Read catchment file ----
# Read one catchment (cat) file and catch possible errors
readcat <- function(f) {
  out <- tryCatch(
    {
    read.table(f, sep="\t", stringsAsFactors = FALSE, header=TRUE,
                 colClasses = c("character", "numeric") ) 
    },
    error=function(cond) {
      message("Error in reading file:", f, "; ", cond, "Continuing without that file.")
      return(NA)
    },
    warning=function(cond) {
      message("Warning in reading file:", f, "; ", cond, "Continuing without that file.")
      return(NA)
    },
    finally={}  # Finally: nothing else.
   )    
  return(out)
}

# Areal average -----

# This function computes the areal average of monthly ETR
#    from the model's cat_*.txt output files.
#    In order to use it as first/final function in mops,
#    it has to use the arguments "parameters, gof, moreArgs" 
process_catfiles <- function(parameters, gof, args) {
  # Find output files and initialize containers
  message(Sys.time(), "(UTC): Post-processing of cat-files with ETR.")
  files = Sys.glob(args["pattern"] )
  files = sort(files)
  vals = c()
  colnames = c()
  # Iterate over files
  for (f in files) {
    ##print(f)
    cat = readcat(f)
    if (!is.data.frame(cat)) {
      next
    }
    if (length(vals)==0) {
      vals = cat[,args["varcol"]]
      dtimes = as.POSIXct(cat[,args["timecol"]], format="%Y-%m-%d %H:%M:%S")
    } else {
      vals = cbind(vals, cat[,args["varcol"]])
    }
    colnames = c(colnames, tools::file_path_sans_ext(f)) 
  }
  names(vals) = colnames
  # Units: from m/s (over an hour) to mm
  vals = vals * 1e3 * 3600
  # Compute average value over all subcatchments
  avg = apply(X=vals, MARGIN=1, FUN=mean, na.rm=TRUE) # mean? oder besser gewichtetes Mittel?
  
  # Compute monthly sums
  months = strftime(dtimes, format="%Y-%m-01 00:00:00")
  monthly = aggregate(avg, by=list(months), "sum", na.rm=TRUE)
  names(monthly) = c("month", "etr")
  
  # Write results file
  ##monthlyfile <- berryFunctions::newFilename(args["outfile"]) # don't overwrite file not possible in second call objfunc
  write.table(monthly, args["outfile"], quote=FALSE, sep = "\t", row.names=FALSE)
  
  gc(verbose=FALSE)
  return(0)
}



# Compute areal average ETR from MODIS data
#   This will produce the "observed" data.frame 
#   to be used with error2 in objfunc.
etr_from_modis <- function(infile="modis_etr/data.txt") {
  cats = read.table(infile, sep="\t", stringsAsFactors = FALSE, header=TRUE,
                    na.strings = "nan" )
  # Compute average over all catchment
  avg = apply(X=cats[,2:ncol(cats)], MARGIN=1, FUN=mean, na.rm=TRUE) # gewichtetes Mittel?
  df = data.frame(month=as.POSIXct(cats$datetime), etr=avg)
  gc(verbose=FALSE)
  return(df)
}


# Objective functions ----

# model_gof, through modelError_multiDim 
# - updates parameter file, 
# - calls the model and extracts the results, 
# - does NOT delete txt, log and html files in output folder outdir (objfunc does that)
# - returns GOF
model_gof <- function(
  my_parameters  = parameters,
  model_path     = mymodelpath,
  my_model_args  = model_args,
  my_func_first  = updatePars_byTemplate,
  moreArgs_first = model_args,
  func_final     = dummygc,
  moreArgs_final = NULL,
  sim_files      = sim_file,
  sim_colsTime   = sim_colTime,
  sim_colsValue  = sim_colValue,
  observed       = obs_file,
  my_obs_colTime = obs_colTime
) {
  # Choose a return value in case of error
  errorout <- list(c(mNSE=9999, pBias=9999, rmse=9999, NSE=-9999))
  names(errorout) <- sim_colValue
  out = tryCatch(
    {
    mops::modelError_multiDim(
      parameters=my_parameters, model_path=model_path, model_args=my_model_args, 
      func_first=my_func_first, moreArgs_first=moreArgs_first, 
      func_final=func_final, moreArgs_final = process_catfiles_args,
      sim_files=sim_files, sim_colsTime=sim_colsTime, sim_colsValue=sim_colsValue,
      observed=observed, obs_colTime=my_obs_colTime,
      sim_colsep     = "\t",
      sim_timeConv   = function(x) {as.POSIXct(strptime(x,"%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "UTC") },
      periods.select = periods,
      periods.ignore = periods.ignore,
      gof_function   = gof_func
      )
    },
    error=function(cond) {
      message("An error occured while calling ECHSE; ", cond)
      return(errorout)
    },
    warning=function(cond) {
      message("An warning occured while calling ECHSE: ", cond)
      return(errorout)
    },
    finally={} # Finally: nothing else.
  )
  return(out)
}




# objective functions considering discharge (and ETR)

objfuncQ <- function(paramsQ, keepmonthly=FALSE) 
            objfunc(paramsQ, gofETR=FALSE, keepmonthly=keepmonthly)

objfunc <- function(params, gofETR=TRUE, keepmonthly=FALSE) {
  begintime <- Sys.time()
  # format parameters as needed by HBV subroutine PARINN
  # sprintf returns a character vector containing a formatted combination of text and variable values
  if(parameter_fmt!="") params = sprintf(parameter_fmt, params)
  
  # Unlog logarithmized parameters for use in hypsoRR -> warum erst hier? 
  params["rate_inter"] = 10.^params["rate_inter"]
  params["rate_base"] = 10.^params["rate_base"]
  params["str_base"] = 10.^params["str_base"]
  
  # set parameter names in parameter vector - this is so dangerous!
  # attr(params, "names")=read.table(file=myrangetable, header=TRUE, sep="\t", colClasses= c("character","numeric","numeric"))[,1] 
  print(params) # control output to console
  message(Sys.time(),"(UTC): objfunc now calls mops::modelError_multiDim for discharge...")

  # 1st call for discharge
  error1 = model_gof(my_parameters  = params,
                     func_final     = process_catfiles,
                     moreArgs_final = process_catfiles_args)
  error1 = round(error1[["qx_avg"]], 5)
  
  gc(verbose=FALSE)
  
  if(gofETR) {
  # 2nd call for ETR
  if(error1["mNSE"]<9999) {
    message(Sys.time(),"(UTC): objfunc now calls mops::modelError_multiDim for ETR...")
    error2 = model_gof(                              # instead of :
      my_parameters  = params,                       #
      model_path     = "calibration/dummymodel.cmd", # mymodelpath "D:/Masterarbeit/Modellierung/Echse/echse_engines/bin/hypsoRR"
      my_model_args  = NULL,                         # model_args (defined at file beginning)
      my_func_first  = dummygc,                      # updatePars_byTemplate
      moreArgs_first = NULL,                         # model_args
      sim_files      = c(etr=as.vector(process_catfiles_args["outfile"])), # sim_file
      sim_colsTime   = c("month"),                   # sim_colTime
      sim_colsValue  = c("etr"),                     # sim_colValue
      observed       = etr_from_modis(),             # obs_file
      my_obs_colTime = "month"                       # obs_colTime
    )
  } else {
    message(Sys.time(),"(UTC): mNSE Discharge was 9999, objfunc now sets ETR gof measures to (-)9999 ...")
    error2 = list(etr=c(mNSE=9999, pBias=9999, rmse=9999, NSE=-9999))
  }
  error2 <-  round(error2[["etr"]], 5)
  } # end if(gofETR)
  
  # final GOF measure + printing
  error <- if(gofETR) total_error(error1, error2) else error1["mNSE"]
  printout <- data.frame(GOF_discharge=error1)
  if(gofETR) printout <- cbind(printout, GOF_ETR=error2)
  print(printout)
  message("total_error = ", error)
  
  # delete files:
  mfile <- process_catfiles_args["outfile"]
  if(keepmonthly) file.copy(mfile, berryFunctions::newFilename(mfile, quiet=TRUE))
  file.remove(mfile) 
  qfile <- process_catfiles_args["outfileq"]
  if(keepmonthly) file.copy(sim_file, berryFunctions::newFilename(qfile, quiet=TRUE))
  empty_outputfolder()
  message(Sys.time(),"(UTC): objfunc is done after ", 
          round(difftime(Sys.time(), begintime, units="min"),2), 
          " minutes. Happyness be to you and your beloved ;-) ...")
  # go home:
  gc(verbose=FALSE)
  return(error)
}


