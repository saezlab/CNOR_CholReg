rm(list=ls())
library(MEIGOR)
library(CNORode)
source("Rscripts/run_opt_pert_functions.R")

seed=as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31);
set.seed(seed);

# read in the data
model <- readSIF(file.path("1_input_data","network.txt"));

cno_data_drugs <-readMIDAS(file.path("./1_input_data","pert_Hela_drugs_prot.csv"))
cno_data_siRNA <-readMIDAS(file.path("./1_input_data","pert_Hela_siRNA_prot.csv"))
cno_data_drugs_metab <-readMIDAS(file.path("./1_input_data","pert_Hela_drugs_metab.csv"))

# merge the data
cno_data_b <- merge.cnodata(cno_data_drugs, cno_data_drugs_metab)
cno_data_c <- merge.cnodata(cno_data_b, cno_data_siRNA)

conditions <- colnames(cno_data_c$dataMatrix[,cno_data_c$TR])
#conditions.rm <- conditions[c(3,4,10:15,20:25)]
#conditions.rm <- conditions[c(3,4,22,23)]
signals <- gsub("DV:", "", colnames(cno_data_c$dataMatrix[,cno_data_c$DV]))
#signals.rm <- signals[!(signals %in% model$namesSpecies)]
#signals.rm <- c(signals.rm, "Q9UBM7", "Q15392") # remove proteins that don't show a difference in levels between LPDS and LPDS + statin


cno_data_c <- simplify.DA(cno_data_c, 50)
#cno_data2_c <- reduce.cnodata(cno_data_c, signals.rm = signals.rm)
#cno_data2_c <- cno_data_c
cno_data_c <- bootstrap_midas_data(cno_data_c)
cnolist <-makePeterCNOlist(cno_data_c)

#templist=cnolist;

cnolist <- reordering.full(cnolist)
cnolist <- addTimePoints(cnolist)

#plotModel(model,cnolist)

ode_parameters <- prepare_ode_parameters(cnolist, model)

  ## Parameter Optimization on Cluster
  #essm
  paramsSSm=defaultParametersSSm()
  paramsSSm$local_solver = "DHC"
  paramsSSm$maxtime = Inf
  paramsSSm$maxeval = 50000;
  paramsSSm$dim_refset=40;
  paramsSSm$ndiverse=500;
  paramsSSm$atol=1e-6;
  paramsSSm$reltol=1e-6;
  paramsSSm$maxStepSize=1;
  paramsSSm$transfer_function = 2;
  paramsSSm$maxNumSteps=10000;
  paramsSSm$nan_fac=1e10;
  
  opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm)
  
  f=opt_pars$ssm_results$f;
  time=opt_pars$ssm_results$time;
  fbest=opt_pars$ssm_results$fbest;
  xbest=opt_pars$ssm_results$xbest;
  neval=opt_pars$ssm_results$neval;
  numeval=opt_pars$ssm_results$numeval;
  
  save(list=c('f','time','fbest','xbest','neval','numeval','seed','model','opt_pars'),file=file.path("./3_results",paste('opt_Hela_',seed,sep='')));
  
  sessionInfo()
