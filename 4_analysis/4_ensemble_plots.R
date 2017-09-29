####
# Script to generate Ensemble plots for the cnodata
####
library(ggplot2)
library(reshape2)

####################################
# Functions
####################################
#' Rename experiments
#' 
#' Renames experiment from E1 to e.g. control for the different data cubes
#' 
#' @param cube.cellline list of data cubes which is output by the function Sim_data_cube
#' 
#' @return cube.cellline data cubes with the renamed experiments
rename_experiments <- function(cube.cellline){
  col.order.exps <- c("control", "02GW", "1GW", "02T09", "1T09", "05HC", "1HC", "2statin", "10statin",
                      "LPDS", "LPDS1statin", "LPDS5statin", 
                      "sSREBP1_s129", "sSREBP1_s130", "sSREBP2_s27", "sSREBP2_s28", "sSREBP1/2", "sLDLR_s06", "sLDLR_s07",
                      "sNPC1_sc", "sNPC1_s69", "sHMGCS1_s62", "sHMGCS1_s63")
  
  for(i in c(1:3)){
    cube.cellline[[i]]$Experiment <- gsub("^E1$", "02GW", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E2$", "02T09", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E3$", "05HC", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E4$", "10statin", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E5$", "1GW", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E6$", "1HC", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E7$", "1T09", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E8$", "2statin", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E9$", "control", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E10$", "LPDS", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E11$", "LPDS1statin", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E12$", "LPDS5statin", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E13$", "sHMGCS1_s62", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E14$", "sHMGCS1_s63", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E15$", "sLDLR_s06", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E16$", "sLDLR_s07", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E17$", "sNPC1_sc", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E18$", "sNPC1_s69", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E19$", "sSREBP1/2", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E20$", "sSREBP1_s129", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E21$", "sSREBP1_s130", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E22$", "sSREBP2_s27", cube.cellline[[i]]$Experiment)
    cube.cellline[[i]]$Experiment <- gsub("^E23$", "sSREBP2_s28", cube.cellline[[i]]$Experiment)
    
    cube.cellline[[i]]$Experiment <- factor(cube.cellline[[i]]$Experiment, levels = col.order.exps)
    
  }
  
  return(cube.cellline)
}

#' Rename Timepoints
#' 
#' Renames timepoints from T1 to e.g. control for the different data cubes
#' 
#' @param cube.cellline list of data cubes which is output by the function Sim_data_cube
#' 
#' @return cube.cellline data cubes with the renamed experiments
rename_timepoints <- function(cube.cellline){
  col.order <- c("T0", "Tx1", "Tx2", "Tx3", "T1")
  
  for(i in c(1:3)){
    cube.cellline[[i]]$Time <- gsub("^T1$", "T0", cube.cellline[[i]]$Time)
    cube.cellline[[i]]$Time <- gsub("^T2$", "Tx1", cube.cellline[[i]]$Time)
    cube.cellline[[i]]$Time <- gsub("^T3$", "Tx2", cube.cellline[[i]]$Time)
    cube.cellline[[i]]$Time <- gsub("^T4$", "Tx3", cube.cellline[[i]]$Time)
    cube.cellline[[i]]$Time <- gsub("^T5$", "T1", cube.cellline[[i]]$Time)
 
    cube.cellline[[i]]$Time <- factor(cube.cellline[[i]]$Time, levels = col.order)
    
  }
  return(cube.cellline)
}


#' Rename nodes
#' 
#' Renames nodes from e.g. Q01581 to HMGCS1
#' 
#' @param cube.cellline list of data cubes which is output by the function Sim_data_cube
#' 
#' @return cube.cellline data cubes with the renamed experiments
rename_nodes <- function(cube.cellline){
  # col.order <- c("control", "02GW", "1GW", "02T09", "1T09", "05HC", "1HC", "2statin", "10statin")

  for(i in c(1:3)){
    cube.cellline[[i]]$Node <- gsub("P14324", "FDPS", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P53396", "ACLY", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q13907", "IDI1", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q15738", "NSDHL", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P48449", "LSS", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P53602", "MVD", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P37268", "FDFT1", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P49327", "FASN", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q13085", "ACACA", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P10515", "DLAT", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("O95573", "ACSL3", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q9BWD1", "ACAT2", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q01581", "HMGCS1", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("P04035", "HMGCR", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q9UBM7", "DHCR7", cube.cellline[[i]]$Node)
    cube.cellline[[i]]$Node <- gsub("Q15392", "DHCR24", cube.cellline[[i]]$Node)
    
    #cube.cellline[[i]]$Node <- factor(cube.cellline[[i]]$Node, levels = col.order)
    
  }
  
  return(cube.cellline)
}


#' Generate data of all the simulation data
#' 
#' It generates the data from a large list of simulation results, the cnolist, and model
#' 
#' @param model_sim_subset large list of model simulation results
#' @param model model from CNORode
#' @param cnolist cnolist object from CNORode
#' 
#' @return cube list of three cubes (predictions, experimental data, perturbations)
#' 
Sim_data_cube <- function(model_sim_subset, model, cnolist){
  # generate cube for the different times
  n.times <- 5
  probs <- c(0, 0.1, 0.25, 0.75, 0.9, 1)
  n.nodes <- length(model$namesSpecies)
  
  # define format of different arrays
  predcube.mean <- array(NA, dim=c(n.times,n.nodes,23), dimnames = list(c(paste("T", seq(1:n.times), sep="")), c(model$namesSpecies), c(paste("E", seq(1:23), sep=""))))
  predcube.exp <- predcube.mean
  predcube.pert <- array(NA, dim=c(n.times,length(cnolist$namesCues),23), dimnames = list(c(paste("T", seq(1:n.times), sep="")), 
                                                                                               c(cnolist$namesStimuli, paste("si", cnolist$namesInhibitors, sep="")), 
                                                                                               c(paste("E", seq(1:23), sep=""))))
  predcube.quant <- array(NA, dim=c(n.times,n.nodes,23, length(probs)), dimnames = list(c(paste("T", seq(1:n.times), sep="")), c(model$namesSpecies), c(paste("E", seq(1:23), sep="")),
                                                                                   c(paste("Q", probs, sep=""))))
  
  # loop over the simulation data and retrieve the data
  for(t in 1:n.times){
    for(n in 1:n.nodes){
      for(e in 1:23){
        mean <- mean(unlist(lapply(model_sim_subset, function(x) x[[t]][e,n])))
        probquants <- quantile(unlist(lapply(model_sim_subset, function(x) x[[t]][e,n])), probs = probs, na.rm=TRUE)
        predcube.mean[t, n, e] <- mean
        predcube.quant[t, n, e, ] <- probquants
      }
    }
  }
  
  # retrieve the experimental data
  for(t in 1:n.times){
    predcube.exp[t, match(cnolist$namesSignals, colnames(predcube.exp)),] <- t(cnolist$valueSignals[[t]])
    predcube.pert[t, ,] <- t(cnolist$valueCues)
  }
  
  # Combine the 
  predtable.mean <- as.data.frame.table(predcube.mean, responseName = "value")
  predtable.quants <- as.data.frame.table(predcube.quant, responseName = "value")
  colnames(predtable.mean) <- c("Time", "Node", "Experiment", "Value")
  colnames(predtable.quants) <- c("Time", "Node", "Experiment", "Quantile", "Value")
  predtable.quants <- dcast(predtable.quants, Time + Node + Experiment ~ Quantile, value.var = "Value")
  
  predtable <- merge(predtable.mean, predtable.quants, by=c("Time", "Node", "Experiment"))
  predtable <- predtable[!is.na(predtable[,4]),]
  predtable <- droplevels(predtable)
  head(predtable)
  
  predtable.exp <- as.data.frame.table(predcube.exp, responseName = "value")
  colnames(predtable.exp) <- c("Time", "Node", "Experiment", "Value.exp")
  predtable.exp <- predtable.exp[!is.na(predtable.exp[,4]),]
  predtable.exp <- droplevels(predtable.exp)
  head(predtable.exp)
  
  predtable.pert <- as.data.frame.table(predcube.pert, responseName = "value")
  colnames(predtable.pert) <- c("Time", "Node", "Experiment", "Value.pert")
  
  # move perturbations nodes to front
  col.order.nodes <- c(levels(predtable.pert$Node), levels(predtable$Node)[!(levels(predtable$Node) %in% levels(predtable.pert$Node))])
  predtable.pert$Node <- factor(predtable.pert$Node, levels = col.order.nodes)
  predtable.exp$Node <- factor(predtable.exp$Node, levels = col.order.nodes)
  
  # only select measured timepoint
  predtable.exp <- subset(predtable.exp, Time == "T5")
  #predtable.sub <- subset(cube.Hela[[1]], Node %in% c("N1", "N5", "N6", "N7")& Experiment %in% c("E1", "E2", "E3", "E4"))
  
  cube.cellline <- list(predtable, predtable.exp, predtable.pert)
  
  cube.cellline <- rename_experiments(cube.cellline)
  cube.cellline <- rename_nodes(cube.cellline)
  #cube.cellline <- rename_timepoints(cube.cellline)
  
  return(cube.cellline)
}

#' Plot an ensemble of functions
#' 
#' It uses the data cube containing the data for a range of ensembl functions to make a plot showing the result of these different fuctions.
#' 
#' @param cube.cellline input data in the format of a list of three tables
#' @return \item{p} ggplot2 plot
#' 
plot_cube <- function(cube.cellline, title){
  p <- ggplot(data = cube.cellline[[1]], aes(x=Time, y=Value, group = Node)) + facet_grid(Experiment ~ Node) 
  p <- p + geom_line() + geom_ribbon(aes(ymin=Q0.25, ymax=Q0.75), alpha=0.5) + geom_ribbon(aes(ymin=Q0.1, ymax=Q0.9), alpha=0.2)
  p <- p + geom_point(data = cube.cellline[[2]], aes(x=Time, y=Value.exp, group=Node), color="red", size=2)
  p <- p + theme(strip.text.y = element_text(angle = 0, hjust = 0, size=15), strip.text.x = element_text(angle = 90, vjust = 0, size=15),
                 axis.title = element_text(size = 15), title = element_text(size = 25), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 legend.title = element_text(size = 12))
  p <- p + labs(title = title) + ylab("Abundance [a.u.]")
  p <- p + geom_rect(data = calculate_error(cube.cellline), aes(group = Node, fill = error),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) 
  p <- p +  scale_fill_gradient2(limits=c(0,1), midpoint = 0.5, low="white", mid = "yellow", high="red", na.value = "grey80")
  p <- p +  scale_y_continuous(breaks=c(0,0.5,1))
  
    # with adding the perturbations
  #p + geom_line(data = cube.cellline[[3]], aes(x=Time, y=Value.pert, group=Node), color="red")
  return(p)
}

plot_cube_selection <- function(cube.celllines, title, celllines, Nodes, Experiments, facets= Experiment ~ Cellline + Node){
  for(i in 1:length(cube.celllines)){
    cube.cellline <- cube.celllines[[i]]
    pred.sel <- subset(cube.cellline[[1]], (Node %in% Nodes) & (Experiment %in% Experiments))
    exp.sel <- subset(cube.cellline[[2]], (Node %in% Nodes) & (Experiment %in% Experiments))
    error.sel <- calculate_error(list(pred.sel, exp.sel))
    
    pred.sel$Cellline <- celllines[i]
    exp.sel$Cellline <- celllines[i]
    error.sel$Cellline <- celllines[i]
    
    if(i == 1){
      pred.sel2 <- pred.sel
      exp.sel2 <- exp.sel
      error.sel2 <- error.sel
    }
    if(i >= 1){
      pred.sel2 <- rbind(pred.sel2, pred.sel)
      exp.sel2 <- rbind(exp.sel2, exp.sel)
      error.sel2 <- rbind(error.sel2, error.sel)
    }
  }
  pred.sel2$Node <- factor(pred.sel2$Node, levels= Nodes)
  exp.sel2$Node <- factor(exp.sel2$Node, levels= Nodes)
  error.sel2$Node <- factor(error.sel2$Node, levels= Nodes)

  pred.sel2$Experiment <- factor(pred.sel2$Experiment, levels= Experiments)
  exp.sel2$Experiment <- factor(exp.sel2$Experiment, levels= Experiments)
  error.sel2$Experiment <- factor(error.sel2$Experiment, levels= Experiments)
  
  pred.sel2$Cellline <- factor(pred.sel2$Cellline, levels= celllines)
  exp.sel2$Cellline <- factor(exp.sel2$Cellline, levels= celllines)
  error.sel2$Cellline <- factor(error.sel2$Cellline, levels= celllines)
  
  p <- ggplot(data = pred.sel2, aes(x=Time, y=Value, group = Node)) + facet_grid(facets) 
  p <- p + geom_line() + geom_ribbon(aes(ymin=Q0.25, ymax=Q0.75), alpha=0.5) + geom_ribbon(aes(ymin=Q0.1, ymax=Q0.9), alpha=0.2)
  p <- p + geom_point(data = exp.sel2, aes(x=Time, y=Value.exp, group=Node), color="red", size=2)
  p <- p + theme(strip.text.y = element_text(angle = 0, hjust = 0, size=15), strip.text.x = element_text(angle = 90, vjust = 0, size=15),
                 axis.title = element_text(size = 15), title = element_text(size = 25), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 legend.title = element_text(size = 12))
  p <- p + labs(title = title) + ylab("Abundance [a.u.]")
  p <- p + geom_rect(data = error.sel2, aes(group = Node, fill = error),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.3) 
  p <- p +  scale_fill_gradient2(limits=c(0,1), midpoint = 0.5, low="white", mid = "yellow", high="red", na.value = "grey80")
  p <- p +  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1))
  
  # with adding the perturbations
  #p + geom_line(data = cube.cellline[[3]], aes(x=Time, y=Value.pert, group=Node), color="red")
  return(p)
}

calculate_error <- function(cube.cellline){
  table.cellline <- merge(cube.cellline[[1]], cube.cellline[[2]], by=c("Time", "Node", "Experiment"), all.x=TRUE)
  table.cellline$error <- abs(table.cellline$Value.exp - table.cellline$Value)
  table.cellline <- table.cellline[,c("Time", "Node", "Experiment", "Value", "Value.exp", "error")]
  table.cellline <- table.cellline[table.cellline$Time == "T5",]
  return(table.cellline)
}

calculate_RMSE <- function(model_sim_subset, model, cnolist){
  # retrieve data from timepoint 5
  exp.data <- cnolist$valueSignals[[5]]
  RMSE <- rep(NA, times = length(model_sim_subset))
  
  for(i in 1:length(model_sim_subset)){
    pred.data <- model_sim_subset[i][[1]][[5]]
    pred.data2 <- pred.data[,match(cnolist$namesSignals, model$namesSpecies)]
    diff.data <- c(pred.data2 - exp.data)
    diff.data <- diff.data[!is.na(diff.data)]
    RMSE[i] <- sqrt(sum((diff.data)^2)/length(diff.data))
  }
  return(RMSE)
}

#############################################################
# Generate data cubes for each cell line with all the results
#############################################################
cube.Hela <- Sim_data_cube(model_sim[ids.Type.Hela], model, cnolist_Hela)
cube.Huh7 <- Sim_data_cube(model_sim[ids.Type.Huh7], model, cnolist_Huh7)
cube.HepG2 <- Sim_data_cube(model_sim[ids.Type.HepG2], model, cnolist_HepG2)
cube.Hek <- Sim_data_cube(model_sim[ids.Type.Hek], model, cnolist_Hek)

# cube when using HepG2 parameters but Huh7 experimental data
cube.HepG2.Huh7 <- Sim_data_cube(model_sim[ids.Type.HepG2], model, cnolist_Huh7)
cube.Huh7.HepG2 <- Sim_data_cube(model_sim[ids.Type.Huh7], model, cnolist_HepG2)

cube.celllines <- list(cube.Hela, cube.Huh7, cube.HepG2, cube.Hek)
celllines = c("Hela", "Huh7", "HepG2", "Hek")

unique(cube.Hela[[1]]$Node)
unique(cube.Hela[[1]]$Experiment)

#example plots
plot_cube_selection(cube.celllines, "FASN vs HMGCS1", celllines = celllines, Nodes = c("FASN", "HMGCS1"), Experiments = c("control", "LPDS1statin", "LPDS5statin"))
plot_cube_selection(list(cube.Huh7, cube.HepG2), "LXR regulation of FASN", celllines = c("Huh7", "HepG2"), Nodes = c("T090137", "GW3965", "FASN", "FDFT1"), Experiments = c("control", "02GW", "1GW", "02T09", "1T09"))
plot_cube_selection(list(cube.Huh7, cube.HepG2), "LXR regulation of FASN", celllines = c("Huh7", "HepG2"), Nodes = c("T090137", "FASN", "FDFT1", "DLAT"), Experiments = c("control", "02T09", "1T09"), facets= Node ~ Cellline + Experiment)
plot_cube_selection(list(cube.Huh7, cube.HepG2), "LXR regulation of FASN", celllines = c("Huh7", "HepG2"), Nodes = c("T090137", "GW3965", "FASN", "FDFT1", "DLAT"), Experiments = c("control", "02GW", "1GW", "02T09", "1T09"), facets= Node ~ Cellline + Experiment)


pdf(paste(Sys.Date(), "Comparison_plots.pdf"),width = 7,height = 7)
plot_cube_selection(list(cube.Huh7, cube.HepG2), "LXR regulation of FASN", celllines = c("Huh7", "HepG2"), Nodes = c("GW3965", "FASN", "FDFT1", "DLAT"), Experiments = c("control", "02GW", "1GW"), facets= Node ~ Cellline + Experiment)
plot_cube_selection(list(cube.Huh7, cube.HepG2), "statin regulation", celllines = c("Huh7", "HepG2"), Nodes = c("atorvastatin", "FASN",  "HMGCS1", "FDFT1"), Experiments = c("control", "2statin", "10statin"), facets= Node ~ Cellline + Experiment)
dev.off()

pdf(paste(Sys.Date(), "Comparison_plots2.pdf"),width = 7,height = 7)
plot_cube_selection(list(cube.Huh7, cube.HepG2), "LXR regulation of FASN", celllines = c("Huh7", "HepG2"), Nodes = c("T090137", "FASN", "FDFT1", "DLAT"), Experiments = c("control", "02T09", "1T09"), facets= Node ~ Cellline + Experiment)
plot_cube_selection(list(cube.Huh7, cube.HepG2), "statin regulation", celllines = c("Huh7", "HepG2"), Nodes = c("atorvastatin", "FASN",  "FDFT1", "HMGCS1"), Experiments = c("control", "2statin", "10statin"), facets= Node ~ Cellline + Experiment)
dev.off()

# Plot Ensembl plots
pdf(paste(Sys.Date(), "Ensemble_plots.pdf"),width = 20,height = 10)
plot_cube(cube.Hela, "Hela")
plot_cube(cube.Huh7, "Huh7")
plot_cube(cube.HepG2, "HepG2")
plot_cube(cube.Hek, "Hek")
plot_cube(cube.HepG2.Huh7, "HepG2 model with Huh7 data")
plot_cube(cube.Huh7.HepG2, "Huh7 model with HepG2 data")
dev.off()


# calculate error
RMSE.Hela <- calculate_RMSE(model_sim[ids.Type.Hela], model, cnolist_Hela)
RMSE.Huh7 <- calculate_RMSE(model_sim[ids.Type.Huh7], model, cnolist_Huh7)
RMSE.HepG2 <- calculate_RMSE(model_sim[ids.Type.HepG2], model, cnolist_HepG2)
RMSE.Hek <- calculate_RMSE(model_sim[ids.Type.Hek], model, cnolist_Hek)

# see how many were above 
RMSE.thr <- 0.1
sum(RMSE.Hela < RMSE.thr)/length(RMSE.Hela)
sum(RMSE.Huh7 < RMSE.thr)/length(RMSE.Huh7)
sum(RMSE.HepG2 < RMSE.thr)/length(RMSE.HepG2)
sum(RMSE.Hek < RMSE.thr)/length(RMSE.Hek)

mean(c(RMSE.Hela, RMSE.Huh7, RMSE.HepG2, RMSE.Hek))
median(c(RMSE.Hela, RMSE.Huh7, RMSE.HepG2, RMSE.Hek))

