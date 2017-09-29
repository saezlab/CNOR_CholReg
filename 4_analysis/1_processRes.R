rm(list=ls())
# set working directory to source file directory. 
seed=as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31);
set.seed(seed);
library(MEIGOR)
library(CNORode)
source("../2_cluster_files/Rscripts/run_opt_pert_functions.R")

#################################################################
# read in results from parameter optimization from results folder
#################################################################
dir.files <- '../3_results'
files=dir(dir.files);

counter=0;
FBEST=c();
XBEST=c();
TYPE=c();

for(i in 1:length(files)){

  if(grepl('opt_Hela', files[i])){
    counter=counter+1;
    load(file.path(dir.files,files[i]));
    FBEST[counter]=fbest;
    XBEST=rbind(XBEST,xbest);
    TYPE=c(TYPE,1);
  }

  if(grepl('opt_HepG2', files[i])){
    counter=counter+1;
    load(file.path(dir.files,files[i]));
    FBEST[counter]=fbest;
    XBEST=rbind(XBEST,xbest);
    TYPE=c(TYPE,2);
  }

  if(grepl('opt_Huh7', files[i])){
    counter=counter+1;
    load(file.path(dir.files,files[i]));
    FBEST[counter]=fbest;
    XBEST=rbind(XBEST,xbest);
    TYPE=c(TYPE,3);
  }

  if(grepl('opt_HEK', files[i])){
    counter=counter+1;
    load(file.path(dir.files,files[i]));
    FBEST[counter]=fbest;
    XBEST=rbind(XBEST,xbest);
    TYPE=c(TYPE,4);
  }
}

# discard results that are above a certain threshold if they exist
FBEST.thr <- 0.05
hist(FBEST, breaks=30)
abline(v = FBEST.thr, col="red", lwd = 3)
sum(FBEST>FBEST.thr)/length(FBEST)
length(FBEST)
TYPE[which(FBEST>FBEST.thr)]

if(sum(FBEST>FBEST.thr) > 0){
  XBEST=XBEST[-which(FBEST>FBEST.thr),];
  TYPE=TYPE[-which(FBEST>FBEST.thr)];
  FBEST=FBEST[-which(FBEST>FBEST.thr)];
}
hist(FBEST)


#############################################################
# read the experimental data and prior-knowledge network
#############################################################

#HELA
#############################################################

# read the data
model <- readSIF(file.path("..","1_input_data","network.txt"));
cno_data_drugs <-readMIDAS(file.path("..","1_input_data","pert_Hela_drugs_prot.csv"))
cno_data_siRNA <-readMIDAS(file.path("..","1_input_data","pert_Hela_siRNA_prot.csv"))
cno_data_drugs_metab <-readMIDAS(file.path("..","1_input_data","pert_Hela_drugs_metab.csv"))

# merge data
cno_data_b <- merge.cnodata(cno_data_drugs, cno_data_drugs_metab)
cno_data_c <- merge.cnodata(cno_data_b, cno_data_siRNA)
cno_data_c <- simplify.DA(cno_data_c, 50)
cnolist <-makePeterCNOlist(cno_data_c)

cnolist <- reordering.full(cnolist)
cnolist_Hela <- addTimePoints(cnolist)


#HEPG2
###############################################################

# read the data
model <- readSIF(file.path("..","1_input_data","network.txt"));
cno_data_drugs <-readMIDAS(file.path("..","1_input_data","pert_HepG2_drugs_prot.csv"))
cno_data_siRNA <-readMIDAS(file.path("..","1_input_data","pert_HepG2_siRNA_prot.csv"))
cno_data_drugs_metab <-readMIDAS(file.path("..","1_input_data","pert_HepG2_drugs_metab.csv"))

# merge data
cno_data_b <- merge.cnodata(cno_data_drugs, cno_data_drugs_metab)
cno_data_c <- merge.cnodata(cno_data_b, cno_data_siRNA)
cno_data_c <- simplify.DA(cno_data_c, 50)
cnolist <-makePeterCNOlist(cno_data_c)

cnolist <- reordering.full(cnolist)
cnolist_HepG2 <- addTimePoints(cnolist)


#Huh7
###############################################################

# read the data
model <- readSIF(file.path("..","1_input_data","network.txt"));
cno_data_drugs <-readMIDAS(file.path("..","1_input_data","pert_Huh7_drugs_prot.csv"))
cno_data_siRNA <-readMIDAS(file.path("..","1_input_data","pert_Huh7_siRNA_prot.csv"))
cno_data_drugs_metab <-readMIDAS(file.path("..","1_input_data","pert_Huh7_drugs_metab.csv"))

# merge data
cno_data_b <- merge.cnodata(cno_data_drugs, cno_data_drugs_metab)
cno_data_c <- merge.cnodata(cno_data_b, cno_data_siRNA)
cno_data_c <- simplify.DA(cno_data_c, 50)
cnolist <-makePeterCNOlist(cno_data_c)

cnolist <- reordering.full(cnolist)
cnolist_Huh7 <- addTimePoints(cnolist)


#HEK
###############################################################
# read the data
model <- readSIF(file.path("..","1_input_data","network.txt"));
cno_data_drugs <-readMIDAS(file.path("..","1_input_data","pert_Hek_drugs_prot.csv"))
cno_data_siRNA <-readMIDAS(file.path("..","1_input_data","pert_Hek_siRNA_prot.csv"))
cno_data_drugs_metab <-readMIDAS(file.path("..","1_input_data","pert_Hek_drugs_metab.csv"))

# merge data
cno_data_b <- merge.cnodata(cno_data_drugs, cno_data_drugs_metab)
cno_data_c <- merge.cnodata(cno_data_b, cno_data_siRNA)
cno_data_c <- simplify.DA(cno_data_c, 50)
cnolist <-makePeterCNOlist(cno_data_c)

cnolist <- reordering.full(cnolist)
cnolist_Hek <- addTimePoints(cnolist)

cnolist <- cnolist_Hela

# combine all cnolist data into one list
cnolist_list <- list(cnolist_Hek, cnolist_Hela, cnolist_HepG2, cnolist_Huh7)

# plot prior-knowledge model
############################
plotModel(model,cnolist_Hek)


# write parameters in one file
#####################################
ode_parameters <- prepare_ode_parameters(cnolist, model)

# obtain the identifiers for each cell line
ids.Type.Hela <- which(TYPE == 1)
ids.Type.HepG2 <- which(TYPE == 2)
ids.Type.Huh7 <- which(TYPE == 3)
ids.Type.Hek <- which(TYPE == 4)

cluster.maxStepSize <- 1
cluster.transfer_fun <- 2
cluster.reltol <- 1e-6
cluster.atol <- 1e-6

model_sim=list();
ode_parameter_list=list();

for(i in 1:nrow(XBEST)){
  #print(i)
  ode_parameter_list[[i]] <- ode_parameters
  ode_parameter_list[[i]]$parValues[ode_parameters$index_opt_pars]=XBEST[i,]
}

##################################
# Simulate how each model behaves
##################################


# calculate the values for the different nodes with the different parameters sets
#################################################################################

print("Hela")
for(i in ids.Type.Hela){
  #print(i);
  model_sim[[i]]=getLBodeModelSim(cnolist_Hela, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
                                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
}

print("HepG2")
for(i in ids.Type.HepG2){
  #print(i);
  model_sim[[i]]=getLBodeModelSim(cnolist_HepG2, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
                                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
}

print("Huh7")
for(i in ids.Type.Huh7){
  #print(i);
  model_sim[[i]]=getLBodeModelSim(cnolist_Huh7, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
                                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
}

print("Hek")
for(i in ids.Type.Hek){
  #print(i);
  ode_parameters$parValues[ode_parameters$index_opt_pars]=XBEST[i,];
  model_sim[[i]]=getLBodeModelSim(cnolist_Hek, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
                                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
}


# print fit for all parameter sets (this takes some time)
#########################################################
# pdf(paste(Sys.Date(), "simulations_all.pdf", sep="_"),width=30,height=15)
# 
# print("Hela")
# plot.new()
# mtext("\nHela",cex=10, adj=0, padj=1)
# for(i in ids.Type.Hela){
#   #print(i);
#   plotLBodeModelSim(cnolist_Hela, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                     maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
#   plotLBodeFitness(cnolist_Hela, model,ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                    maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
# }
# 
# print("HepG2")
# plot.new()
# mtext("\nHepG2",cex=10, adj=0, padj=1)
# for(i in ids.Type.HepG2){
#   #print(i);
#   plotLBodeModelSim(cnolist_HepG2, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                     maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
#   plotLBodeFitness(cnolist_HepG2, model,ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                    maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
# }
# 
# print("Huh7")
# plot.new()
# mtext("\nHuh7",cex=10, adj=0, padj=1)
# for(i in ids.Type.Huh7){
#   #print(i);
#   plotLBodeModelSim(cnolist_Huh7, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                     maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
#   plotLBodeFitness(cnolist_Huh7, model,ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                    maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
# }
# 
# print("Hek")
# plot.new()
# mtext("\nHEK",cex=10, adj=0, padj=1)
# for(i in ids.Type.Hek){
#   #print(i);
#   plotLBodeModelSim(cnolist_Hek, model, ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                     maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
#   plotLBodeFitness(cnolist_Hek, model,ode_parameter_list[[i]], reltol = cluster.reltol, atol = cluster.atol, 
#                    maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)
# }
# 
# dev.off()

# print fit for best parameter set
#########################################################

print("simulation")
colnames(XBEST)=ode_parameters$parNames[ode_parameters$index_opt_pars];
xopt=XBEST[which(FBEST==min(FBEST)),];
ode_parameters$parValues[ode_parameters$index_opt_pars]=xbest;


pdf(paste(Sys.Date(), "simulation.pdf", sep="_"),width=30,height=15)
plotModel(model,cnolist_Hela)

  temp=FBEST[TYPE==1];
  best=which(FBEST==min(temp))
  plot.new()
  mtext("\nHela",cex=10, adj=0, padj=1)
  plotLBodeModelSim(cnolist_Hela, model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                    maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
  plotLBodeFitness(cnolist_Hela, model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                   maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)

temp=FBEST[TYPE==2];
best=which(FBEST==min(temp))
  plot.new()
  mtext("\nHepG2",cex=10, adj=0, padj=1)
  plotLBodeModelSim(cnolist_HepG2, model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                    maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
  plotLBodeFitness(cnolist_HepG2, model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                   maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


temp=FBEST[TYPE==3];
best=which(FBEST==min(temp))
  plot.new()
  mtext("\nHuh7",cex=10, adj=0, padj=1)

  plotLBodeModelSim(cnolist_Huh7, model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                    maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
  plotLBodeFitness(cnolist_Huh7, model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                   maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


temp=FBEST[TYPE==4];
best=which(FBEST==min(temp))
  plot.new()
  mtext("\nHEK",cex=10, adj=0, padj=1)

  plotLBodeModelSim(cnolist_Hek, model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                    maxStepSize =cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
  plotLBodeFitness(cnolist_Hek, model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                   maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


dev.off()

# print fit for best parameter set only for two timepoints
#########################################################

pdf(paste(Sys.Date(), "simulation_2timepoints.pdf", sep="_"),width=30,height=15)
plotModel(model,cnolist_Hela)

temp=FBEST[TYPE==1];
best=which(FBEST==min(temp))
plot.new()
mtext("\nHela",cex=10, adj=0, padj=1)
plotLBodeModelSim(rmTimePoints(cnolist_Hela), model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                  maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
plotLBodeFitness(rmTimePoints(cnolist_Hela), model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)

temp=FBEST[TYPE==2];
best=which(FBEST==min(temp))
plot.new()
mtext("\nHepG2",cex=10, adj=0, padj=1)
plotLBodeModelSim(rmTimePoints(cnolist_HepG2), model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                  maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
plotLBodeFitness(rmTimePoints(cnolist_HepG2), model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


temp=FBEST[TYPE==3];
best=which(FBEST==min(temp))
plot.new()
mtext("\nHuh7",cex=10, adj=0, padj=1)

plotLBodeModelSim(rmTimePoints(cnolist_Huh7), model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                  maxStepSize = cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
plotLBodeFitness(rmTimePoints(cnolist_Huh7), model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


temp=FBEST[TYPE==4];
best=which(FBEST==min(temp))
plot.new()
mtext("\nHEK",cex=10, adj=0, padj=1)

plotLBodeModelSim(rmTimePoints(cnolist_Hek), model, ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                  maxStepSize =cluster.maxStepSize,timeSignals=seq(0,50), transfer_function = cluster.transfer_fun)
plotLBodeFitness(rmTimePoints(cnolist_Hek), model,ode_parameter_list[[best]], reltol = cluster.reltol, atol = cluster.atol, 
                 maxStepSize = cluster.maxStepSize, transfer_function = cluster.transfer_fun)


dev.off()

# plot the different parameters
####################################
pdf(paste(Sys.Date(), "pars.pdf", sep="_"))
df.XBEST <- data.frame(XBEST, row.names=c(1:dim(XBEST)[1]))
df.XBEST <- cbind(df.XBEST, data.frame(TYPE = TYPE))
df.XBEST$TYPE <- gsub("^1$", "Hela", df.XBEST$TYPE)
df.XBEST$TYPE <- gsub("^2$", "HepG2", df.XBEST$TYPE)
df.XBEST$TYPE <- gsub("^3$", "Huh7", df.XBEST$TYPE)
df.XBEST$TYPE <- gsub("^4$", "Hek", df.XBEST$TYPE)


#library(vegan)
cv <- apply(df.XBEST[,c(1:dim(df.XBEST)[2]-1)], 2, function(x) sd(x)/mean(x))
head(cv)
cv[order(cv, decreasing=TRUE)]

stdev <- apply(df.XBEST[,c(1:dim(df.XBEST)[2]-1)], 2, sd)
head(stdev)
stdev[order(stdev, decreasing=TRUE)]



# insert parameters that have been fixed
parNames_nonoptimized <- which(!(ode_parameters$parNames %in% colnames(df.XBEST)))
for(i in parNames_nonoptimized){
  df.XBEST[,ode_parameters$parNames[i]] <- ode_parameters$parValues[i]
}

# plot all parameters within the boundaries that were set
for(i in 1:length(ode_parameters$parNames[ode_parameters$index_opt_pars])){
  boundaries <- c(ode_parameters$LB[ode_parameters$index_opt_pars][i], ode_parameters$UB[ode_parameters$index_opt_pars][i])
  p <- ggplot(df.XBEST, aes_string(x="TYPE", y=ode_parameters$parNames[ode_parameters$index_opt_pars][i], group="TYPE")) + geom_boxplot() + ggtitle(ode_parameters$parNames[ode_parameters$index_opt_pars][i])
  p <- p + scale_y_continuous(trans = "log", limits = boundaries)
  print(p)
}

dev.off()

# write data into file 
save.image(file.path(getwd(), paste(Sys.Date(),"results.RData")))
writeLines(capture.output(sessionInfo()), paste(Sys.Date(), "sessionInfo.txt"))