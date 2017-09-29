reordering.full <- function(cnolist){
  if(identical(cnolist$namesStimuli,c("Input05","statin_1uM","statin_2uM","statin_5uM","statin_10uM","LDPS",
                                      "05HC","1HC","02T09","1T09","02GW","1GW"))){
    #### reordering
    
    # "Input05"     "statin_1uM"  "statin_2uM"  "statin_5uM"  "statin_10uM" "LDPS"        "05HC"        "1HC"        
    # "02T09"       "1T09"        "02GW"        "1GW"            
    
    # "statin_1uM"  "statin_2uM"  "statin_5uM"  "statin_10uM
    cnolist$namesStimuli[2]="statin";
    colnames(cnolist$valueStimuli)[2]="statin";
    cnolist$valueStimuli[,2]=cnolist$valueStimuli[,2]*0.3+cnolist$valueStimuli[,3]*0.45+cnolist$valueStimuli[,4]*0.6+cnolist$valueStimuli[,5]*0.9;
    cnolist$valueStimuli=cnolist$valueStimuli[,-c(3,4,5)];
    cnolist$namesStimuli=cnolist$namesStimuli[-c(3,4,5)];
    
    # "Input05" "statin"  "LDPS"    "05HC"    "1HC"     "02T09"   "1T09"    "02GW"    "1GW"    
    
    #"05HC","1HC"
    cnolist$namesStimuli[4]="HC";
    colnames(cnolist$valueStimuli)[4]="HC";
    cnolist$valueStimuli[,4]=cnolist$valueStimuli[,4]*0.4+cnolist$valueStimuli[,5]*0.9;
    cnolist$valueStimuli=cnolist$valueStimuli[,-5];
    cnolist$namesStimuli=cnolist$namesStimuli[-5];
    
    # "Input05" "statin"  "LDPS"    "HC"      "02T09"   "1T09"    "02GW"    "1GW"  
    
    #"02T09","1T09"
    cnolist$namesStimuli[5]="T09";
    colnames(cnolist$valueStimuli)[5]="T09";
    cnolist$valueStimuli[,5]=cnolist$valueStimuli[,5]*0.25+cnolist$valueStimuli[,6]*0.9;
    cnolist$valueStimuli=cnolist$valueStimuli[,-6];
    cnolist$namesStimuli=cnolist$namesStimuli[-6];
    
    # "Input05" "statin"  "LDPS"    "HC"      "T09"     "02GW"    "1GW"    
    
    # "02GW"     "1GW"
    cnolist$namesStimuli[6]="GW";
    colnames(cnolist$valueStimuli)[6]="GW";
    cnolist$valueStimuli[,6]=cnolist$valueStimuli[,6]*0.25+cnolist$valueStimuli[,7]*0.9;
    cnolist$valueStimuli=cnolist$valueStimuli[,-7];
    cnolist$namesStimuli=cnolist$namesStimuli[-7];
  }
  
  # reorder inhibitors
  if(isTRUE(identical(cnolist$namesInhibitors, c("SREBP1_s129","SREBP1_s130","SREBP2_s27","SREBP2_s28","LDLR_s06","LDLR_s07",
                                                 "NPC1_s69","NPC1_c","Q01581_s62","Q01581_s63")))){
    #"SREBP1_s129" "SREBP1_s130" "SREBP2_s27"  "SREBP2_s28"  "LDLR_s06"    "LDLR_s07"    "NPC1_s69"    "NPC1_c"     
    #"Q01581_s62"  "Q01581_s63"    
    cnolist$namesInhibitors=c("SREBP1","SREBP2", "LDLR","NPC1", "Q01581")
    cnolist$valueInhibitors=cbind(apply(cnolist$valueInhibitors[,1:2], 1, max),
                                  apply(cnolist$valueInhibitors[,3:4], 1, max),
                                  apply(cnolist$valueInhibitors[,5:6], 1, max),
                                  apply(cnolist$valueInhibitors[,7:8], 1, max),
                                  apply(cnolist$valueInhibitors[,9:10], 1, max))
    colnames(cnolist$valueInhibitors)=c("SREBP1","SREBP2","LDLR","NPC1", "Q01581")
    
    cnolist$valueInhibitors=cnolist$valueInhibitors*0.6;
    
    cnolist$valueInhibitors[cnolist$valueInhibitors>1]=1;
  }
  
  cnolist$valueStimuli[cnolist$valueStimuli==0]=0.05;
  cnolist$valueStimuli[cnolist$valueStimuli==1]=0.9;
  
  cnolist$namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);
  cnolist$valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
  
  return(cnolist)
  
}



simplify.DA <- function(cno_data, time){
  for(i in cno_data$DAcol){
    cno_data$dataMatrix[which(cno_data$dataMatrix[,i] > 0 ),i] <- time
  }
  return(cno_data)
}

# simplify.DA <- function(cno_data){
#   time.points <- cno_data$dataMatrix[,cno_data$DAcol]
#   time.points <- rowMeans(time.points)
#   time.points[time.points > 0] <- 50
#   colnames(cno_data$dataMatrix)[cno_data$DAcol[1]] <- "DA:ALL"
#   cno_data$dataMatrix[,cno_data$DAcol[1]] <- time.points
#   cno_data$dataMatrix <- cno_data$dataMatrix[,c(1:cno_data$DAcol[1],cno_data$DVcol)]
#   cno_data$DVcol <- cno_data$DVcol - (cno_data$DVcol[1] - cno_data$DAcol[1]) + 1
#   cno_data$DAcol <- c(cno_data$DAcol[1])
#   return(cno_data)
# }


reduce.cnodata <- function(cno_data, conditions.rm=NULL, signals.rm = NULL){
  if(!is.null(conditions.rm)){
    conditions.rm2 <- conditions.rm[!conditions.rm %in% c("TR:InputNPC1", "TR:InputCitrate")] # remove conditions as otherwise all rows are removed
    if(length(conditions.rm2)>0){
      TR.rm <- cno_data$dataMatrix[,conditions.rm2]
      rows.rm <- unique(which(rowSums(TR.rm) > 0))
    }
    TR.cols.rm <- which(colnames(cno_data$dataMatrix) %in% conditions.rm)
    
    if(length(conditions.rm2)>0){
      cno_data$dataMatrix <- cno_data$dataMatrix[-rows.rm,]
    }
    cno_data$dataMatrix <- cno_data$dataMatrix[,-which(colnames(cno_data$dataMatrix) %in% conditions.rm)]

    cno_data$TRcol <- cno_data$TRcol[1:(length(cno_data$TRcol) - length(TR.cols.rm))]
    cno_data$DAcol <- cno_data$DAcol - length(TR.cols.rm)
    cno_data$DVcol <- cno_data$DVcol - length(TR.cols.rm)
  }
  if(!is.null(signals.rm)){
    DV.cols.rm <- which(colnames(cno_data$dataMatrix) %in% paste("DV:",signals.rm, sep=""))
    DA.cols.rm <- which(colnames(cno_data$dataMatrix) %in% paste("DA:",signals.rm, sep=""))
    
    cno_data$dataMatrix <- cno_data$dataMatrix[,-c(DA.cols.rm, DV.cols.rm)]
    cno_data$DAcol <- cno_data$DAcol[1:(length(cno_data$DAcol) - length(DA.cols.rm))]
    cno_data$DVcol <- cno_data$DVcol[1:(length(cno_data$DVcol) - length(DV.cols.rm))] - length(DA.cols.rm)
  }
  return(cno_data)
}

merge.cnodata <- function(cno_data, cno_data2){
  #check columns not yet in data
  columns.added <- colnames(cno_data2$dataMatrix)[!(colnames(cno_data2$dataMatrix) %in% colnames(cno_data$dataMatrix))]
  
  grep("^TR:", columns.added)
  grep("^DA:", columns.added)
  grep("^DV:", columns.added)
  
  TR.columns.added <- columns.added[grep("^TR:", columns.added)]
  DA.columns.added <- columns.added[grep("^DA:", columns.added)]
  DV.columns.added <- columns.added[grep("^DV:", columns.added)]
  
  # number of columns added
  TR.columns.added.n <- length(TR.columns.added)
  DA.columns.added.n <- length(DA.columns.added)
  DV.columns.added.n <- length(DV.columns.added)
  
  # add columns to dataMatrix
  cno_data$dataMatrix[,c(TR.columns.added, DA.columns.added, DV.columns.added)] <- 0
  # indices of data
  TRcol.added <- which(colnames(cno_data$dataMatrix) %in% TR.columns.added)
  DAcol.added <- which(colnames(cno_data$dataMatrix) %in% DA.columns.added)
  DVcol.added <- which(colnames(cno_data$dataMatrix) %in% DV.columns.added)
  
  # merge data into new data.matrix
  #   cno_data$dataMatrix <- cbind(cno_data$dataMatrix[,c(cno_data$TRcol)],
  #                                cno_data$dataMatrix[,TRcol.added],
  #                                cno_data$dataMatrix[,c(cno_data$DAcol)],
  #                                cno_data$dataMatrix[,DAcol.added],
  #                                cno_data$dataMatrix[,c(cno_data$DVcol)],
  #                                cno_data$dataMatrix[,DVcol.added])
  #   head(cno_data$dataMatrix)
  
  # merge DA and DV - for adding metabolites
  if(identical(cno_data$dataMatrix[,cno_data$TRcol], cno_data2$dataMatrix[,cno_data2$TRcol]) & length(TRcol.added) == 0){
    cno_data$dataMatrix <- cbind(cno_data$dataMatrix[,c(cno_data$TRcol)],
                                 cno_data$dataMatrix[,c(cno_data$DAcol)],
                                 cno_data2$dataMatrix[,DA.columns.added],
                                 cno_data$dataMatrix[,c(cno_data$DVcol)],
                                 cno_data2$dataMatrix[,DV.columns.added])
  }
  
  # merge TR - for adding the siRNA data
  if(length(TRcol.added) > 0 & length(DAcol.added) == 0 & length(DVcol.added) == 0 ){
    cno_data$dataMatrix <- cbind(cno_data$dataMatrix[,c(cno_data$TRcol)],
                                 cno_data$dataMatrix[,TRcol.added],
                                 cno_data$dataMatrix[,c(cno_data$DAcol)],
                                 cno_data$dataMatrix[,c(cno_data$DVcol)])
    
    # add additional data.
    # add columns to imported data.frame to allow rbind
    columns.added.rev <- colnames(cno_data$dataMatrix)[!(colnames(cno_data$dataMatrix) %in% colnames(cno_data2$dataMatrix))]
    DA.columns.added.rev <- columns.added.rev[grep("^DA:", columns.added.rev)]
    DV.columns.added.rev <- columns.added.rev[grep("^DV:", columns.added.rev)]
    DA.columns.added.rev.n <- length(DA.columns.added.rev)
    
    cno_data2$dataMatrix[,c(DA.columns.added.rev)] <- data.frame(do.call(cbind, rep(list(cno_data2$dataMatrix[,cno_data2$DAcol[1]]), DA.columns.added.rev.n)))
    cno_data2$dataMatrix[,c(DV.columns.added.rev)] <- NA
    
    cno_data2$dataMatrix <- cno_data2$dataMatrix[,colnames(cno_data$dataMatrix)]
    
    cno_data$dataMatrix <- rbind(cno_data$dataMatrix, cno_data2$dataMatrix)
    
  }
  
  
  TRcol.added2 <- which(colnames(cno_data$dataMatrix) %in% TR.columns.added)
  DAcol.added2 <- which(colnames(cno_data$dataMatrix) %in% DA.columns.added)
  DVcol.added2 <- which(colnames(cno_data$dataMatrix) %in% DV.columns.added)
  
  cno_data$TRcol <- c(cno_data$TRcol, TRcol.added2)
  cno_data$DAcol <- c(cno_data$DAcol + TR.columns.added.n, DAcol.added2)
  cno_data$DVcol <- c(cno_data$DVcol + TR.columns.added.n + DA.columns.added.n, DVcol.added2)
  
  return(cno_data)
}

fix.parameters <- function(ode_parameters, edges, value){
  ode_parameters$parValues[which(ode_parameters$parNames %in% edges)] <- value
  k <- which(ode_parameters$parNames %in% edges)
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% k)]
  return(ode_parameters)
}



prepare_ode_parameters <- function(cnolist, model){
  ## set initial parameters
  ode_parameters=createLBodeContPars(model, LB_n = 0.5, LB_k = 0.1, LB_tau = 0.5,
                                     UB_n = 4, UB_k = 10, UB_tau = 5, default_n = 2,
                                     default_k = 0.54, default_tau = 2, opt_n = TRUE, opt_k = TRUE,
                                     opt_tau = TRUE, random = FALSE)
  
  
  modelSim=plotLBodeModelSim(cnolist, model, ode_parameters, timeSignals=seq(0,2,0.2), verbose=TRUE, transfer_function = 2)
  
   
  ### mass-action rate constants
  ###########################################

  print("WARNING: THIS ORDER VARIES IF THE SIF VARIES. Uncomment conde in rhsODE.c to check new indices!!!");
  
  print(paste("k1 is equivalent to:", ode_parameters$parNames[44]))
  print(paste("k2 is equivalent to:", ode_parameters$parNames[49]))
  print(paste("k3 is equivalent to:", ode_parameters$parNames[60]))
  print(paste("k4 is equivalent to:", ode_parameters$parNames[69]))
  print(paste("k5 is equivalent to:", ode_parameters$parNames[13]))
  print(paste("k6 is equivalent to:", ode_parameters$parNames[7]))
  print(paste("k7 is equivalent to:", ode_parameters$parNames[47]))
  print(paste("k8 is equivalent to:", ode_parameters$parNames[9]))

  print(paste("kM4 is equivalent to:", ode_parameters$parNames[68]))

  print(paste("kIatorvastatin is equivalent to:", ode_parameters$parNames[65]))

  print(paste("tauAcetyl_CoA_index is equivalent to:", ode_parameters$parNames[62]))
  print(paste("tauAcetoacetyl_CoA_index is equivalent to:", ode_parameters$parNames[48]))
  print(paste("tauHMG_CoA_index is equivalent to:", ode_parameters$parNames[58]))
  print(paste("tauMevalonate_index is equivalent to:", ode_parameters$parNames[64]))
  print(paste("tauCholER_index is equivalent to:", ode_parameters$parNames[6]))
  
  
  # collect the parameters that are used for mass-action kinetics
  met.react <- c("Input05_k_Acetyl_CoA", "Acetyl_CoA_k_Acetoacetyl_CoA", "Acetoacetyl_CoA_k_HMG_CoA", 
                 "HMG_CoA_k_Mevalonate", "Mevalonate_k_CholER", 
                 "CholMedia_k_CholER", "ACAT2act_k_Acetoacetyl_CoA", "LDLR_k_CholER")
  
  kM.react <- c("HMG_CoA_n_Mevalonate")
  
  kI.react <- c("atorvastatin_k_Mevalonate")
  
  tau.react <- c("HMGCS1act_k_HMG_CoA", "Acetyl_CoA_n_Acetoacetyl_CoA", "Acetyl_CoA_k_HMG_CoA", 
                  "atorvastatin_n_Mevalonate", "CholMedia_n_CholER")
  
  tau.react.fast <- c("CholMedia_n_CholER")
 

  #adjust bounds for the metabolic reactions
  ode_parameters$UB[which(ode_parameters$parNames %in% met.react)] <- 0.1
  ode_parameters$LB[which(ode_parameters$parNames %in% met.react)] <- 1e-3
  ode_parameters$parValues[which(ode_parameters$parNames %in% met.react)] <- 0.01
  
  met.react.Acetyl_CoA <- c("Input05_k_Acetyl_CoA", "ACAT2act_k_Acetoacetyl_CoA")

  ode_parameters$UB[which(ode_parameters$parNames %in% met.react.Acetyl_CoA)] <- 1e6
  ode_parameters$LB[which(ode_parameters$parNames %in% met.react.Acetyl_CoA)] <- 1e5
  ode_parameters$parValues[which(ode_parameters$parNames %in% met.react.Acetyl_CoA)] <- 1e5
  
  ode_parameters$LB[which(ode_parameters$parNames %in% c("Acetyl_CoA_k_Acetoacetyl_CoA"))] <- 1e-5
  
  #adjust bounds for kM
  ode_parameters$parValues[which(ode_parameters$parNames == "HMG_CoA_n_Mevalonate")] <- 0.5
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% which(ode_parameters$parNames %in% c("HMG_CoA_n_Mevalonate")))]
  
  #adjust bounds for kI
  ode_parameters$UB[which(ode_parameters$parNames %in% kI.react)] <- 0.01
  ode_parameters$LB[which(ode_parameters$parNames %in% kI.react)] <- 0.0001
  ode_parameters$parValues[which(ode_parameters$parNames %in% kI.react)] <- 0.0003

  #tau parameters for metabolic reactions are fixed, only allowing the k to vary
  ode_parameters$UB[which(ode_parameters$parNames %in% tau.react)] <- 10
  ode_parameters$LB[which(ode_parameters$parNames %in% tau.react)] <- 0.1
  ode_parameters$parValues[which(ode_parameters$parNames %in% tau.react)] <- 1
  
  ode_parameters$UB[which(ode_parameters$parNames %in% tau.react.fast)] <- 300
  ode_parameters$LB[which(ode_parameters$parNames %in% tau.react.fast)] <- 100
  ode_parameters$parValues[which(ode_parameters$parNames %in% tau.react.fast)] <- 100
    
  k <- which(ode_parameters$parNames %in% tau.react)
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% k)]
  

  ### drug-uptake and drug effect
  ###########################################

  # drug uptake making it a sigmoidal curve n > 2
  ode_parameters$LB[which(ode_parameters$parNames == "T09_n_T090137")] <- 2
  ode_parameters$LB[which(ode_parameters$parNames == "GW_n_GW3965")] <- 2
  ode_parameters$LB[which(ode_parameters$parNames == "statin_n_atorvastatin")] <- 2


  # adjust bounds to ensure minimal effect on intracellular drug concentration and transcription factor
  ode_parameters$LB[which(ode_parameters$parNames == "tau_GW3965")] <- 1
  ode_parameters$LB[which(ode_parameters$parNames == "tau_T090137")] <- 1
  ode_parameters$LB[which(ode_parameters$parNames == "tau_atorvastatin")] <- 1

  ode_parameters$UB[which(ode_parameters$parNames == "T090137_k_LXR")] <- 0.5
  ode_parameters$UB[which(ode_parameters$parNames == "GW3965_k_LXR")]  <- 0.5
  ode_parameters$UB[which(ode_parameters$parNames == "T090137_n_LXR")] <- 2
  ode_parameters$UB[which(ode_parameters$parNames == "GW3965_n_LXR")]  <- 2

  
  # fix tau for these reactions to 1
  ode_parameters <- fix.parameters(ode_parameters, c("tau_atorvastatin", "tau_T090137", "tau_GW3965", "tau_LXR"), 1)

  # fix HC_n_SREBP edge to 2
  ode_parameters <- fix.parameters(ode_parameters, c('HC_n_SREBP'), 2)

  # allow drug metabolites a wider range to approximate a higher background value
  ode_parameters$LB[which(ode_parameters$parNames %in% c("atorvastatin_n_Atorvastatin_Lactone", "atorvastatin_n_hydroxy_atorvastatin", "GW3965_n_hydroxy_GW3965"))] <- 0.1
  ode_parameters$UB[which(ode_parameters$parNames %in% c("atorvastatin_k_Atorvastatin_Lactone", "atorvastatin_k_hydroxy_atorvastatin", "GW3965_k_hydroxy_GW3965"))] <- 100


  # fix LDPS on CholMedia
  ode_parameters <- fix.parameters(ode_parameters, c('LDPS_n_CholMedia'), 1)
  ode_parameters <- fix.parameters(ode_parameters,c('LDPS_k_CholMedia'), 0.2)
  ode_parameters <- fix.parameters(ode_parameters, c('Input05_n_CholMedia'), 3)
  ode_parameters <- fix.parameters(ode_parameters, c('Input05_k_CholMedia'), 0.5)
  ode_parameters <- fix.parameters(ode_parameters, c('tau_CholMedia'), 1)

  
  ### logic-based intracellular interactions
  ###########################################

  #adjust bounds for fast reactions such as SREBP response and activity of enzymes
  ode_parameters$UB[which(ode_parameters$parNames %in% c("tau_SREBP"))] <- 1000
  ode_parameters$LB[which(ode_parameters$parNames %in% c("tau_SREBP"))] <- 100
  ode_parameters$parValues[which(ode_parameters$parNames %in% c("tau_SREBP"))] <- 800
   
  ode_parameters$UB[which(ode_parameters$parNames %in% c("tau_HMGCS1act", "tau_ACAT2act"))] <- 100
  ode_parameters$LB[which(ode_parameters$parNames %in% c("tau_HMGCS1act", "tau_ACAT2act"))] <- 20
  ode_parameters$parValues[which(ode_parameters$parNames %in% c("tau_HMGCS1act", "tau_ACAT2act"))] <- 50

 
 
  # fix SREBP1 and SREBP2 edges from SREBP
  ode_parameters <- fix.parameters(ode_parameters, c("SREBP_n_SREBP2", "SREBP_n_SREBP1"), 3)
  ode_parameters <- fix.parameters(ode_parameters, c("SREBP_k_SREBP2", "SREBP_k_SREBP1"), 0.5)
  ode_parameters <- fix.parameters(ode_parameters, c("tau_SREBP2", "tau_SREBP1"), 800)

  # fix CholER_n_SREBP edge according to Radhakrishnan publication to 3.7
   ode_parameters <- fix.parameters(ode_parameters, c('CholER_n_SREBP'), 3.7)
 
  # adjust bounds for nodes to ensure that they are active
  ode_parameters$UB[which(ode_parameters$parNames == "SREBP2_k_P04035")] <- 1
  ode_parameters$UB[which(ode_parameters$parNames == "SREBP2_k_LDLR")] <- 1
  ode_parameters$UB[which(ode_parameters$parNames == "Q9BWD1_k_ACAT2act")] <- 1
  ode_parameters$UB[which(ode_parameters$parNames == "Q01581_k_HMGCS1act")] <- 1
  ode_parameters$LB[which(ode_parameters$parNames == "CholER_n_SREBP")] <- 3
  ode_parameters$LB[which(ode_parameters$parNames == "tau_LDLR")] <- 2



  # fix edges to remove correlations between parameters
  ########################################

  # fix Input05_n_NPC1 edge to 1.5 to only allow variation in Input05_k_NPC1
  ode_parameters <- fix.parameters(ode_parameters, c('Input05_n_NPC1'), 1.5)
  
  # fix SREBP2_n_LDLR to allow variation only in SREBP2_k_LDLR 
  ode_parameters <- fix.parameters(ode_parameters, c('SREBP2_n_LDLR'), 3.5)

  # fix tau parameter for HMGCR
  ode_parameters <- fix.parameters(ode_parameters, c("tau_P04035"), 5)

  # remove edges from optimization that are not used in network
  ########################################

  k <- which(ode_parameters$parNames %in% c("tau_CholER", "tau_Acetoacetyl_CoA", "tau_Acetyl_CoA", 'tau_HMG_CoA', 'tau_Mevalonate'))
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% k)]
  k <- which(ode_parameters$parNames %in% c("Input05_n_Acetyl_CoA", "Mevalonate_n_Acetyl_CoA"))
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% k)]
  k <- which(ode_parameters$parNames %in% c("P04035_n_Mevalonate", "P04035_k_Mevalonate"))
  ode_parameters$index_opt_pars <- ode_parameters$index_opt_pars[-which(ode_parameters$index_opt_pars %in% k)]
  
  return(ode_parameters)
}


makePeterCNOlist<-function(cno_data2_c){
  
  cnolist=list();
  
  cno_data2_c$dataMatrix=as.matrix(cno_data2_c$dataMatrix, rownames.force = TRUE);
  
  cno_data2_c$dataMatrix[is.nan(cno_data2_c$dataMatrix)]=NA;
  
  identifier_table=c();
  for(i in 1:dim(cno_data2_c$dataMatrix)[1]){
    identifier_table=c(identifier_table,paste(cno_data2_c$dataMatrix[i,cno_data2_c$TRcol],sep='',collapse = ''))
  }
  
  unique_identifiers=unique(identifier_table);
  timeSignals=sort(unique(cno_data2_c$dataMatrix[,cno_data2_c$DAcol][,1]));
  
  valueSignals=list();
  valueVariances=list();
  
  for(i in 1:length(timeSignals)){
    valueSignals[[i]]=matrix(NA,length(unique_identifiers),length(cno_data2_c$DVcol));
    valueVariances[[i]]=matrix(0.1,length(unique_identifiers),length(cno_data2_c$DVcol));                      
  }
  
  treatment_names=colnames(cno_data2_c$dataMatrix[,cno_data2_c$TRcol]);
  inhibitorCols=grep(':inhibitor',treatment_names);
  stimuliCols=1:length(treatment_names);
  if(length(inhibitorCols > 0)){
    stimuliCols=stimuliCols[-inhibitorCols];
  }
  namesStimuli=treatment_names[stimuliCols];
  namesInhibitors=treatment_names[inhibitorCols];
  
  valueInhibitors=matrix(NA,length(unique_identifiers),length(inhibitorCols));
  valueStimuli=matrix(NA,length(unique_identifiers),length(stimuliCols));
  
  for(i in 1:dim(cno_data2_c$dataMatrix)[1]){
    rowVals=cno_data2_c$dataMatrix[i,];
    stimuliVals=rowVals[cno_data2_c$TRcol[stimuliCols]];
    inhibitorVals=rowVals[cno_data2_c$TRcol[inhibitorCols]];
    signalVals=rowVals[cno_data2_c$DVcol];
    rowN=which(unique_identifiers==identifier_table[i]);
    timeN=which(timeSignals==rowVals[cno_data2_c$DAcol[1]]);
    valueSignals[[timeN]][rowN,]=colMeans(rbind(valueSignals[[timeN]][rowN,],signalVals),na.rm = TRUE);
    valueStimuli[rowN,]=colMeans(rbind(valueStimuli[rowN,],stimuliVals),na.rm = TRUE);
    valueInhibitors[rowN,]=colMeans(rbind(valueInhibitors[rowN,],inhibitorVals),na.rm = TRUE);
    valueSignals[[timeN]][is.nan(valueSignals[[timeN]])]=NA;
  }
  
  for(i in 1:length(unique_identifiers)){
    rowVals=cno_data2_c$dataMatrix[identifier_table==unique_identifiers[i], , drop=FALSE];
    rowVals=as.matrix(rowVals);
    
    signalVals=rowVals[,cno_data2_c$DVcol];
    if(length(dim(signalVals))<2){
      signalVals=as.matrix(t(signalVals));
      
    }
    
    for(j in 1:length(timeSignals)){
      nSignals=dim(valueSignals[[j]])[2];
      for(k in 1:nSignals){
        #print(signalVals)
        timeN=which(rowVals[,cno_data2_c$DAcol[1]]==timeSignals[j]);
        #print(timeN)
        if(length(timeN)>1){
          #print(signalVals[timeN,k])
          valueVariances[[j]][i,k]=var(signalVals[timeN,k],na.rm = TRUE)
        }else{
          valueVariances[[j]][i,k]=0;
        }
      }
    }
    
    
  }
  
  for(i in 1:dim(valueSignals[[1]])[1]){
    valueSignals[[1]][i,]=colMeans(valueSignals[[1]],na.rm = TRUE);
    
  }
  
  cnolist$valueSignals=valueSignals;
  cnolist$valueStimuli=valueStimuli;
  cnolist$valueInhibitors=valueInhibitors;
  cnolist$valueVariances=valueVariances;
  cnolist$timeSignals=timeSignals;
  cnolist$namesSignals=colnames(cno_data2_c$dataMatrix[,cno_data2_c$DVcol])
  cnolist$namesSignals=sub("DV:",'', cnolist$namesSignals);
  cnolist$namesStimuli=colnames(cno_data2_c$dataMatrix[,cno_data2_c$TRcol])[stimuliCols];
  cnolist$namesStimuli=sub("TR:",'', cnolist$namesStimuli);
  cnolist$namesInhibitors=colnames(cno_data2_c$dataMatrix[,cno_data2_c$TRcol])[inhibitorCols];
  cnolist$namesInhibitors=sub("TR:",'',cnolist$namesInhibitors);
  cnolist$namesInhibitors=sub(":inhibitor",'',cnolist$namesInhibitors);
  cnolist$namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);
  cnolist$valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
  
  colnames(cnolist$valueStimuli)=cnolist$namesStimuli;
  colnames(cnolist$valueInhibitors)=cnolist$namesInhibitors;
  
  
  return(cnolist);
  
}

bootstrap_midas_data<-function(cno_data2_c){
  
  identifier_table=c();
  for(i in 1:dim(cno_data2_c$dataMatrix)[1]){
    identifier_table=c(identifier_table,paste(cno_data2_c$dataMatrix[i,cno_data2_c$TRcol],cno_data2_c$dataMatrix[i,cno_data2_c$DAcol],sep='',collapse = ''))
  }
  
  
  unique_identifiers=unique(identifier_table);
  mat=matrix(NA,length(unique_identifiers),length(cno_data2_c$DVcol));
  
  for(i in 1:length(cno_data2_c$DVcol)){
    lines=which(!is.na(cno_data2_c$dataMatrix[,cno_data2_c$DVcol[i]]));
    lines=sample(lines,replace = TRUE)
    
    vals=list();
    for(j in 1:length(unique_identifiers)){
      vals[[j]]=c(NA);
    }
    
    for (j in 1:length(lines)){
      index_unique=NA;
      for (k in 1:length(unique_identifiers)){
        
        if(identical(unique_identifiers[k],identifier_table[lines[j]])) {
          index_unique=k;
        }
      }
      
      if(is.na(index_unique))stop("something is wrong");
      #print(index_unique)
      #print(lines[j])
      vals[[index_unique]]=c(vals[[index_unique]],cno_data2_c$dataMatrix[lines[j],cno_data2_c$DVcol[i]]);    
    }
    for(j in 1:length(unique_identifiers)){
      mat[j,i]=mean(vals[[j]],na.rm=TRUE);
    }
  }
  mat2=c();
  for(i in 1:length(unique_identifiers)){
    index=which(identifier_table==unique_identifiers[i]);
    mat2=rbind(mat2,c(cno_data2_c$dataMatrix[index[1],cno_data2_c$TRcol],cno_data2_c$dataMatrix[index[1],cno_data2_c$DAcol]))
  }
  mat2 <- apply(mat2, 2, as.numeric)
  
  mat[is.nan(mat)]=NA;
  #colnames(mat)=colnames(cno_data2_c$dataMatrix[,cno_data2_c$DAcol])
  colnames(mat)=colnames(cno_data2_c$dataMatrix[,cno_data2_c$DVcol])
  #cno_data2_c$dataMatrix=cbind(mat2,mat);
  cno_data2_c$dataMatrix=data.frame(cbind(mat2,mat))
  colnames(cno_data2_c$dataMatrix) <- gsub("([[:alpha:]]{2})\\.(.*)", "\\1:\\2", colnames(cno_data2_c$dataMatrix))
  colnames(cno_data2_c$dataMatrix) <- gsub("\\.inhibitor", ":inhibitor", colnames(cno_data2_c$dataMatrix))
  return(cno_data2_c);
  
}

addTimePoints <- function(cnolist){
  timeSignals <- cnolist$timeSignals
  cnolist$timeSignals <- seq(min(timeSignals), max(timeSignals), length.out = 5)
  
  diffSignals <- cnolist$valueSignals[[2]]-cnolist$valueSignals[[1]]
  cnolist$valueSignals[[5]] <- cnolist$valueSignals[[2]]
  cnolist$valueSignals[[2]] <- cnolist$valueSignals[[1]]+0.5*diffSignals
  cnolist$valueSignals[[3]] <- cnolist$valueSignals[[1]]+0.8*diffSignals
  cnolist$valueSignals[[4]] <- cnolist$valueSignals[[1]]+0.9*diffSignals
  
  cnolist$valueVariances[[3]] <- cnolist$valueVariances[[2]]
  cnolist$valueVariances[[4]] <- cnolist$valueVariances[[2]]
  cnolist$valueVariances[[5]] <- cnolist$valueVariances[[2]]
  
  return(cnolist)

}

rmTimePoints <- function(cnolist){
  timeSignals <- cnolist$timeSignals
  timeIds <- c(1, length(timeSignals))
  cnolist$timeSignals <- cnolist$timeSignals[timeIds]
  
  cnolist$valueSignals[[2]] <- cnolist$valueSignals[[5]]
  cnolist$valueVariances[[2]] <- cnolist$valueVariances[[5]]
  cnolist$valueSignals[[3]]
  cnolist$valueSignals[[4]]
  cnolist$valueSignals[[5]]
  cnolist$valueVariances[[3]]
  cnolist$valueVariances[[4]]
  cnolist$valueVariances[[5]]
  
  return(cnolist)
  
}
