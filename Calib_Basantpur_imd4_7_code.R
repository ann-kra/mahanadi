
# mops_ppse analysis Script Basantpur 2017-09-19

rm(list=ls())

# Working directory and functions ----

setwd("D:/Masterarbeit/Modellierung/Calibration_Ann")
source("calibration/Kalib_Basantpur_imd4_7_funs.R")


# SETTINGS ----

periods        <- data.frame(begin=ISOdatetime(2001,3,2,8,0,0),end=ISOdatetime(2010,12,31,23,0,0))
periods.ignore <- data.frame(begin=ISOdatetime(1988,1,1,0,0,0),end=ISOdatetime(1989,12,31,0,0,0))
outdir         <- "calibration/out"
mymodelpath    <- "D:/Masterarbeit/Modellierung/Echse/echse_engines/bin/hypsoRR" # ggf. ändern in ../Echse/
myrangetable   <- "calibration/tbl_ranges_bigger.txt"
sim_file       <- c(qx_avg="calibration/out/gage_Basantpur.txt")
sim_colTime    <- c("end_of_interval")
sim_colValue   <- c("qx_avg")
obs_file       <- read.table("calibration/flow_Basantpur_0.txt", sep="\t", stringsAsFactors=FALSE, colClasses=c("POSIXct","numeric"), header=TRUE)
obs_colTime    <- c("end_of_interval")
mops_log       <- "calibration/mcs_and_pso.log" 
ppso_log       <- "calibration/dds.log" 
ppso_proj      <- "calibration/dds.pro"  
parameter_fmt  <-  "" # set to "" for ECHSE and "%13.4f" for HBV Nordic
if_fails_clear <-  c("out", "phiMod.err.html", "phiMod.log") #for ECHSE: c("out", "phiMod.err.html", "phiMod.log")

process_catfiles_args <-  c( # Argument vector for process_catfiles
  pattern = "calibration/out/cat_*.txt", 
  timecol = "end_of_interval", 
  varcol  = "etr",
  outfile = "calibration/outMonthly/etr_monthly.txt",
  outfileq= "calibration/outMonthly/q_monthly.txt"
)
model_args     <-  c(
  file_control                   ="calibration/cnf_imd4_Basantpur_Ann.txt",
  file_err                       = paste(outdir,"cnf_raingages.err.html",sep="/"),
  file_log                       = paste(outdir,"cnf_raingages.log",sep="/"),
  format_err                     = "html", 
  silent                         = "true",
  table_objectDeclaration        = "out_csa100_Arc_Gis_splitting/Basantpur/objDecl.gage_Basantpur.txt",
  cat_hypsoRR_numParamsIndividual= "calibration/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
  cat_hypsoRR_numParamsShared    = "../echse_projekt_Ann/data/params/cat_num_shared.txt",
  outputDirectory                = outdir,
  file_template                  = "calibration/paramNum_cat_MCS_template_maik.txt", # vorher paramNum_cat_MCS_template2_new.txt
  file_result                    = "calibration/paramNum_cat_MCS_template_updated_mult_str_surf.txt",
  char_open                      = "{",
  char_close                     = "}"
)

# read parameter ranges
param_bounds <-  read.table(file=myrangetable, 
                            header=TRUE, 
                            sep="\t", 
                            colClasses= c("character","numeric","numeric"))[,2:3]
rownames(param_bounds) <-  read.table(file=myrangetable, 
                                      header=TRUE, 
                                      sep="\t", 
                                      colClasses= c("character","numeric","numeric"))[,1]


# Try things: 
# readcat("calibration/out/cat_460.txt")
# process_catfiles(NULL, NULL, process_catfiles_args)
# plot(as.POSIXct(out$month), out$etr, type="l")


if(FALSE){ #don't run this part when sourcing, because it takes ages to run

# Calibration ----

# Call the optimisation routine from ppso package, using our objective function

  
# ############################
# Berry says: optim_dds tries to *minimize* mNSE / total_error. 
# Maybe objfunc should be changed to return - mNSE ?
# ############################

result = ppso::optim_dds(
  objective_function        = objfuncQ, # objfuncQ to calibrate only on discharge, objfunc to calibrate on both discharge and ETR
  number_of_parameters      = length(param_bounds[,1]),
  number_of_particles       = 1,
  max_number_function_calls = 500,
  r                         = 0.2,
  abstol                    = -Inf, 
  reltol                    = -Inf, 
  max_wait_iterations       = 50,
  parameter_bounds          = param_bounds,
  lhc_init                  = TRUE, 
  do_plot                   = NULL, 
  wait_for_keystroke        = FALSE,
  logfile                   = ppso_log,  
  projectfile               = ppso_proj,
  load_projectfile          = "no", 
  break_file                = NULL, 
  plot_progress             = FALSE, 
  tryCall                   = FALSE)


# plot progress ----
ppso::plot_optimization_progress(
         logfile                = ppso_log,
         projectfile            = ppso_proj,
         progress_plot_filename = NULL,
         goodness_plot_filename = NULL,
         cutoff_quantile        = 0.95,
         verbose                = FALSE)





# produce data for plotting, run objfunc with best parameters ----

best_params <- t(read.table(ppso_proj, header=TRUE, sep="\t")[,1:11])[,1]
names(best_params) <- gsub("best_", "", names(best_params))

#run_with_best_params  <-  objfunc(best_params, keepmonthly=TRUE)
run_with_best_paramsQ <- objfuncQ(best_params, keepmonthly=TRUE)

obs_file$end_of_interval <- as.Date(obs_file$end_of_interval, format="%Y-%m-%d %H:%M:%S")
obs_file_subset <- subset(obs_file)
mfiles <- dir(dirname(process_catfiles_args["outfileq"]), full.names=TRUE); mfiles 
monthly <- read.table(tail(mfiles,1), header=TRUE, stringsAsFactors=FALSE, sep="\t")
monthly$end_of_interval <- as.Date(monthly$end_of_interval, format="%Y-%m-%d %H:%M:%S")

# insert here aggregate - Fun = mean
merged_data <- merge(monthly, obs_file, by="end_of_interval", all.x=T)
merged_data$end_of_interval <- as.Date(merged_data$end_of_interval, format="%Y-%m-%d")
plot(merged_data$end_of_interval, merged_data$qx_avg.y, xaxt="n", type="l", las=1)
lines(merged_data$end_of_interval, merged_data$qx_avg.x, col="red")
berryFunctions::monthAxis(ym=TRUE)

max(obs_file$qx_avg) # BLACK (observed streamflow)
max(monthly$qx_avg) # RED  -> (calibrated streamflow)
} 

