library(gridExtra)
library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

#####################
### Functions
#####################

# Function to plot all the different solutions for a given edge and cell line
plot.mult.functions <- function(source.node, target.node, cellline, df.XBEST, FBEST, inhibiting = FALSE, cnolist=NULL){
  
  # define the starting node and end node
  n.edge <- paste(source.node, "n", target.node, sep="_")
  k.edge <- paste(source.node, "k", target.node, sep="_")
  
  # indices for best and good (best 10%) of estimations
  temp=FBEST[df.XBEST$TYPE==cellline];
  good=which(df.XBEST$TYPE==cellline & FBEST <= quantile(temp, 0.1, names=FALSE)) 
  best=which(FBEST==min(temp))
  
  # n and k parameters
  parameters <- data.frame(n = df.XBEST[df.XBEST$TYPE == cellline, n.edge], k = df.XBEST[df.XBEST$TYPE == cellline, k.edge])
  parameters.fgood <- data.frame(n=df.XBEST[good, n.edge], k = df.XBEST[good, k.edge])
  parameters.fbest <- data.frame(n=df.XBEST[best, n.edge], k = df.XBEST[best, k.edge])
  
  # indices for for loops
  k1 <- nrow(parameters)
  k2 <- nrow(parameters.fgood)
  k3 <- nrow(parameters.fbest)
  
  # create different functions for all, good and best function
  if(inhibiting == FALSE){
    function.list <- list()
    for(i in 1:k1){
      function.list[[i]] <- stat_function(fun=function(x, n, k) {(x^n)/(k^n + x^n)}, 
                                          args=list(n=parameters[i,"n"], k=parameters[i, "k"]), 
                                          colour="grey", size=0.2)
    }
    for(i in (k1+1):(k1+k2)){
      i.k <- i - k1
      function.list[[i]] <- stat_function(fun=function(x, n, k) {(x^n)/(k^n + x^n)}, 
                                          args=list(n=parameters.fgood[i.k,"n"], k=parameters.fgood[i.k, "k"]), 
                                          colour="black", size=0.2)
    }
    for(i in (k1+k2+1):(k1+k2+k3)){
      i.k <- i - (k1 + k2)
      function.list[[i]] <- stat_function(fun=function(x, n, k) {(x^n)/(k^n + x^n)}, 
                                          args=list(n=parameters.fbest[i.k,"n"], k=parameters.fbest[i.k, "k"]), 
                                          colour="red", size=1)
    }
    
  }
  if(inhibiting == TRUE){
    function.list <- list()
    for(i in 1:k1){
      function.list[[i]] <- stat_function(fun=function(x, n, k) {1-((x^n)/(k^n + x^n))}, 
                                          args=list(n=parameters[i,"n"], k=parameters[i, "k"]), 
                                          colour="grey", size=0.2)
    }
    for(i in (k1+1):(k1+k2)){
      i.k <- i - k1
      function.list[[i]] <- stat_function(fun=function(x, n, k) {1-((x^n)/(k^n + x^n))}, 
                                          args=list(n=parameters.fgood[i.k,"n"], k=parameters.fgood[i.k, "k"]), 
                                          colour="black", size=0.2)
    }
    for(i in (k1+k2+1):(k1+k2+k3)){
      i.k <- i - (k1 + k2)
      function.list[[i]] <- stat_function(fun=function(x, n, k) {1-((x^n)/(k^n + x^n))}, 
                                          args=list(n=parameters.fbest[i.k,"n"], k=parameters.fbest[i.k, "k"]), 
                                          colour="red", size=1)
    }
    
  }
  
  # retrieve data if possible
  if(!is.null(cnolist)){
    # find out if source is stimuli or signal and if there is a signal for the target node
    source.stimuli <- sum(source.node %in% cnolist$namesStimuli)
    source.signal <- sum(source.node %in% cnolist$namesSignals)
    target.signal <- sum(target.node %in% cnolist$namesSignals)
    
    # define color depending if source and target has been measured
    if(sum(source.stimuli, source.signal) & target.signal){
      col <- "green4"
    }
    else(col <- "blue4")
    
    # retrieve data from cnolist
    if(source.stimuli){
      l <- which(cnolist$namesStimuli == source.node)
      source.data <- cnolist$valueStimuli[,l]
      source.data.var <- rep(0, length(source.data))
    }
    if(source.signal){
      l <- which(cnolist$namesSignals == source.node)
      source.data <- cnolist$valueSignals[[length(cnolist$valueSignals)]][,l]
      source.data.var <- cnolist$valueVariances[[length(cnolist$valueSignals)]][,l]
    }
    # if there is no data write 0
    if(sum(source.stimuli, source.signal) == 0){
      source.data <- rep(-0, dim(cnolist$valueSignals[[length(cnolist$valueSignals)]])[1])
      source.data.var <- rep(0, length(source.data))
    }
    
    if(target.signal){
      l <- which(cnolist$namesSignals == target.node)
      target.data <- cnolist$valueSignals[[length(cnolist$valueSignals)]][,l]
      target.data.var <- cnolist$valueVariances[[length(cnolist$valueSignals)]][,l]
    }
    # if there is no data write 0
        if(!target.signal){
      target.data <- rep(-0, dim(cnolist$valueSignals[[length(cnolist$valueSignals)]])[1])
      target.data.var <- rep(0, length(target.data))
    }
    
    exp.data <- data.frame(source.data = source.data, target.data = target.data, source.var = source.data.var, target.var = target.data.var)
    exp.data <- exp.data[!is.na(exp.data$source.data),]
    exp.data <- exp.data[!is.na(exp.data$target.data),]
    
  }
  
  
  # plot function
  if(is.null(cnolist)){
    p <- ggplot(data.frame(x=c(0,1)), aes(x)) + function.list + ggtitle(cellline) + xlim(c(0,1)) + ylim(c(0,1))+ xlab(source.node) + ylab(target.node)
    
    # return plot and parameters
    return(list(p, parameters, parameters.fgood, parameters.fbest, exp.data <- NULL))
  }
  
  if(!is.null(cnolist)){
    p <- (ggplot(data.frame(x=c(0,1)), aes(x)) + function.list + ggtitle(cellline) + xlim(c(0,1)) + ylim(c(0,1))
          + xlab(source.node) + ylab(target.node) 
          + geom_point(aes(x=source.data, y=target.data), data = exp.data, size=2, color=col) 
          + geom_errorbar(data = exp.data, aes(ymin=target.data-target.var, 
                                               ymax = target.data+target.var, x=source.data), 
                          width=0, size=1, color=col)
          + geom_errorbarh(data = exp.data, aes(xmin=source.data-source.var, 
                                                xmax = source.data+source.var, x=source.data, y=target.data), 
                           height=0, size=1, color=col)
    )
    
    # return plot and parameters
    return(list(p, parameters, parameters.fgood, parameters.fbest, exp.data))
    }
}

# Function to plot all four cell lines in one plot
plot.4CL.functions <- function(source.node, target.node, df.XBEST, FBEST, inhibiting, cnolist_list){
  plot.HEK <- plot.mult.functions(source.node, target.node, "Hek", df.XBEST, FBEST, inhibiting, cnolist_list[[1]])[[1]]
  plot.Hela <- plot.mult.functions(source.node, target.node, "Hela", df.XBEST, FBEST,inhibiting, cnolist_list[[2]])[[1]]
  plot.HepG2 <- plot.mult.functions(source.node, target.node, "HepG2", df.XBEST, FBEST, inhibiting, cnolist_list[[3]])[[1]]
  plot.Huh7 <- plot.mult.functions(source.node, target.node, "Huh7", df.XBEST, FBEST, inhibiting, cnolist_list[[4]])[[1]]
  
  p <- arrangeGrob(plot.HEK, plot.Hela, plot.HepG2, plot.Huh7, ncol=2,top=paste(source.node, target.node, sep=" > "))
  return(p)
}

# retrieve information if edge is inhibiting
define.inh.edges <- function(edges.Hill.parNames, model){
  edges.Hill.parNames.2 <- strsplit(edges.Hill.parNames, "-")
  edges.Hill.parNames.df  <-  as.data.frame(matrix(unlist(edges.Hill.parNames.2), ncol=2, byrow=T))
  colnames(edges.Hill.parNames.df) <- c("source.node", "target.node")
  edges.Hill.parNames.df <- edges.Hill.parNames.df[order(edges.Hill.parNames.df$source.node, edges.Hill.parNames.df$target.node),]
  edges.Hill.parNames.df$inhibiting <- NA
  
  model.split <- strsplit(model$reacID, split="=")
  for(i in 1:nrow(edges.Hill.parNames.df)){
    
    node.source <- edges.Hill.parNames.df[i,"source.node"]
    node.target <- edges.Hill.parNames.df[i,"target.node"]
    
    l <- grep(node.source, model.split)
    for(k in l){
      m.source <- grep(node.source, model.split[[k]])
      if(sum(m.source %in% 1)){
        m.target <- grep(node.target, model.split[[k]])
        if(sum(m.target %in% 2)){
          m.inh <- grep(paste("!", node.source, sep=""), model.split[[k]])
          if(length(m.inh) > 0){
            edges.Hill.parNames.df[i, "inhibiting"] <- TRUE
          }
          else(edges.Hill.parNames.df[i, "inhibiting"] <- FALSE)
        }
      }
    }
  }
  return(edges.Hill.parNames.df)
}


#####################
# Definitions
#####################
# collect the edges that are used for mass-action kinetics
met.react <- c("Input05_k_Acetyl_CoA", "Acetyl_CoA_k_Acetoacetyl_CoA", "Acetoacetyl_CoA_k_HMG_CoA", 
               "HMG_CoA_k_Mevalonate", "Mevalonate_k_CholER", 
               "CholMedia_k_CholER", "ACAT2act_k_Acetoacetyl_CoA", "LDLR_k_CholER")

kM.react <- c("HMG_CoA_n_Mevalonate")

kI.react <- c("atorvastatin_k_Mevalonate")

NA.react <- c("Input05_n_Acetyl_CoA", "Acetoacetyl_CoA_n_HMG_CoA", "Mevalonate_n_Acetyl_CoA", "P04035_n_Mevalonate", 
              "P04035_k_Mevalonate", "HMGCS1act_k_HMG_CoA", "Acetyl_CoA_n_Acetoacetyl_CoA", "Acetyl_CoA_k_HMG_CoA", 
              "atorvastatin_n_Mevalonate", "CholMedia_n_CholER")

# make dataframe with edges that are described by Hill-type function
edges.Hill <- ode_parameters$parNames[!(ode_parameters$parNames %in% c(met.react, kM.react, kI.react, NA.react))]
edges.Hill <- edges.Hill[edges.Hill %in% ode_parameters$parNames[ode_parameters$index_opt_pars]]
edges.Hill <- edges.Hill[c(grep("_n_", edges.Hill), grep("_k_", edges.Hill))]
edges.Hill <- strsplit(edges.Hill, split = "_")
edges.Hill.n <- lapply(edges.Hill, function(x) {l <- which(x == "n"); if(length(l) > 0){return(paste(paste(x[seq(1,(l-1))], collapse="_"), paste(x[seq((l+1),length(x))], collapse="_"), sep="-"))}})
edges.Hill.k <- lapply(edges.Hill, function(x) {l <- which(x == "k"); if(length(l) > 0){return(paste(paste(x[seq(1,(l-1))], collapse="_"), paste(x[seq((l+1),length(x))], collapse="_"), sep="-"))}})
edges.Hill.parNames <- intersect(unlist(edges.Hill.n), unlist(edges.Hill.k))

setdiff(intersect(unlist(edges.Hill.n), unlist(edges.Hill.k)), edges.Hill.parNames)

edges.Hill.parNames.df <- define.inh.edges(edges.Hill.parNames, model)
edges.Hill.parNames.df

# make dataframe with edges that are described by Hill-type function without removing fixed edges
edges.Hill <- ode_parameters$parNames[!(ode_parameters$parNames %in% c(met.react, kM.react, kI.react, NA.react))]
edges.Hill <- edges.Hill[c(grep("_n_", edges.Hill), grep("_k_", edges.Hill))]
edges.Hill <- strsplit(edges.Hill, split = "_")
edges.Hill.n <- lapply(edges.Hill, function(x) {l <- which(x == "n"); if(length(l) > 0){return(paste(paste(x[seq(1,(l-1))], collapse="_"), paste(x[seq((l+1),length(x))], collapse="_"), sep="-"))}})
edges.Hill.k <- lapply(edges.Hill, function(x) {l <- which(x == "k"); if(length(l) > 0){return(paste(paste(x[seq(1,(l-1))], collapse="_"), paste(x[seq((l+1),length(x))], collapse="_"), sep="-"))}})
edges.Hill.parNames <- intersect(unlist(edges.Hill.n), unlist(edges.Hill.k))

setdiff(intersect(unlist(edges.Hill.n), unlist(edges.Hill.k)), edges.Hill.parNames)

edges.Hill.parNames.df2 <- define.inh.edges(edges.Hill.parNames, model)
edges.Hill.parNames.df2


#####################################################
# Plotting of transfer functions describing the edges
#####################################################

# plot a Hill-type function with certain n and k values
ggplot(data.frame(x=c(0,1)), aes(x)) + stat_function(fun=function(x, n, k) {(x^n)/(k^n + x^n)}, args=list(n=3.99, k=0.45)) + xlim(c(0,1)) + ylim(c(0,1))
ggplot(data.frame(x=c(0,1)), aes(x)) + stat_function(fun=function(x, n, k) {(x^n)/(k^n + x^n)}, args=list(n=1.208, k=1.737)) + xlim(c(0,1)) + ylim(c(0,1))

# plot a certain edge
p <- plot.4CL.functions("statin", "atorvastatin", df.XBEST, FBEST, inhibiting = FALSE, cnolist_list)
grid.draw(p)
ggsave("statin_atorvastatin.pdf", p, width=100, height=70, scale=2, units="mm", dpi=600)

i <- 1
p <- plot.mult.functions(edges.Hill.parNames.df[i, "source.node"], edges.Hill.parNames.df[i, "target.node"], "Huh7", df.XBEST, FBEST, inhibiting = edges.Hill.parNames.df[i, "inhibiting"], cnolist = cnolist_list[[3]])

# plot all transfer functions (time-intensive)
pdf(paste(Sys.Date(), "transfer_functions.pdf"))
for(i in 1:nrow(edges.Hill.parNames.df)){
  #print(i)
  p <- plot.4CL.functions(edges.Hill.parNames.df[i, "source.node"], edges.Hill.parNames.df[i, "target.node"], df.XBEST, FBEST, inhibiting = edges.Hill.parNames.df[i, "inhibiting"], cnolist_list = cnolist_list)
  grid.newpage()
  grid.draw(p)
}
dev.off()

# plot all transfer functions (including fixed ones)
pdf(paste(Sys.Date(), "transfer_functions2.pdf"))
for(i in 1:nrow(edges.Hill.parNames.df2)){
  #print(i)
  p <- plot.4CL.functions(edges.Hill.parNames.df2[i, "source.node"], edges.Hill.parNames.df2[i, "target.node"], df.XBEST, FBEST, inhibiting = edges.Hill.parNames.df2[i, "inhibiting"], cnolist_list = cnolist_list)
  grid.newpage()
  grid.draw(p)
}
dev.off()

#####################################################
# Calculating and plot flux
#####################################################
flux.list <- list()
for(id.condition in 1:dim(cnolist_list[[1]]$valueSignals[[1]])[1]){
  flux.df <- data.frame(Condition = id.condition, Cellline = rep(NA, nrow(df.XBEST)))
  # annotate cell lines
  flux.df[ids.Type.Hek,"Cellline"] <- "HEK"
  flux.df[ids.Type.Hela,"Cellline"] <- "Hela"
  flux.df[ids.Type.Huh7,"Cellline"] <- "Huh7"
  flux.df[ids.Type.HepG2,"Cellline"] <- "HepG2"
  
  flux.df$k1 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "Input05_k_Acetyl_CoA")]
  flux.df$k2 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "Acetyl_CoA_k_Acetoacetyl_CoA")]
  flux.df$k3 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "Acetoacetyl_CoA_k_HMG_CoA")]
  flux.df$k4 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "HMG_CoA_k_Mevalonate")]
  flux.df$k5 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "Mevalonate_k_CholER")]
  flux.df$k6 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "CholMedia_k_CholER")]
  flux.df$k7 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "ACAT2act_k_Acetoacetyl_CoA")]
  flux.df$k8 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "LDLR_k_CholER")]
  
  #flux.df$kM4 <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "HMG_CoA_n_Mevalonate")]
  flux.df$kIatorvastatin <- df.XBEST[,which(ode_parameters$parNames[ode_parameters$index_opt_pars] == "atorvastatin_k_Mevalonate")]
  
  # Flux Input -> Acetyl_CoA
  flux.df$Input_AcetylCoA <- flux.df$k1 * 0.5
  
  # Flux Acetyl_CoA -> Acetoacetyl_CoA
  flux.df$ACAT2act <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "ACAT2act")]))
  flux.df$AcetylCoA <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "Acetyl_CoA")]))
  flux.df$AcetylCoA_AcetoacetylCoA <- flux.df$k2 * flux.df$ACAT2act *flux.df$AcetylCoA^2
  
  # Flux AcetoacetylCoA -> HMG_CoA
  flux.df$HMGCS1act <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "HMGCS1act")]))
  flux.df$AcetoacetylCoA <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "Acetoacetyl_CoA")]))
  flux.df$AcetoacetylCoA_HMGCoA <- flux.df$k3 * flux.df$HMGCS1act *flux.df$AcetylCoA* flux.df$AcetoacetylCoA
  
  # Flux HMG_CoA -> Mevalonate
  flux.df$HMGCR <-  unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "P04035")]))
  flux.df$HMGCoA <-  unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "HMG_CoA")]))
  flux.df$atorvastatin <-  unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "atorvastatin")]))
  # flux.df$HMGCoA_Mev <- flux.df$k4 * flux.df$HMGCR *flux.df$HMGCoA / (flux.df$kM4 + flux.df$HMGCoA + ((flux.df$kM4*flux.df$atorvastatin)/flux.df$kIatorvastatin))
  # flux.df$HMGCRInh <- (flux.df$kM4*flux.df$atorvastatin)/flux.df$kIatorvastatin
  flux.df$HMGCoA_Mev <- flux.df$k4 * flux.df$HMGCR *flux.df$HMGCoA / (0.5 + flux.df$HMGCoA + ((0.5*flux.df$atorvastatin)/flux.df$kIatorvastatin))
  flux.df$HMGCRInh <- (0.5*flux.df$atorvastatin)/flux.df$kIatorvastatin
  
  # Flux Mevalonate -> CholER
  flux.df$Mev = unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "Mevalonate")]))
  flux.df$Mev_Chol <- flux.df$k5 * flux.df$Mev
  
  # Flux Cholesterol uptake
  flux.df$CholMedia <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "CholMedia")]))
  flux.df$LDLR <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "LDLR")]))
  flux.df$NPC1 <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "NPC1")]))
  flux.df$Choluptake <- flux.df$k6 * flux.df$CholMedia * flux.df$LDLR * flux.df$NPC1
  
  # Flux Cholesterol shunt
  flux.df$CholER <- unlist(lapply(model_sim, function(x) x[[5]][id.condition,which(model$namesSpecies == "CholER")]))
  flux.df$Cholshunt <- flux.df$k8 * flux.df$CholER
  
  # Flux Acetyl_CoA shunt
  flux.df$AcetylCoAshunt <- flux.df$k7 * flux.df$AcetylCoA
  
  flux.list[[id.condition]] <- flux.df
}

flux.table <- ldply(flux.list, data.frame)
flux.table$Condition <- factor(flux.table$Condition)

head(flux.df)

# corresponds to cno_list row number
id.condition <- 9 #untreated

# plot fluxes
pdf(paste(Sys.Date(), "fluxes.pdf"))
ggplot(flux.list[[id.condition]], aes(y=Input_AcetylCoA, x=Cellline)) + geom_boxplot() + ggtitle("Input -> Acetyl-CoA")
ggplot(flux.list[[id.condition]], aes(y=AcetylCoA_AcetoacetylCoA, x=Cellline)) + geom_boxplot() + ggtitle("Acetyl-CoA -> Acetoacetyl-CoA")
ggplot(flux.list[[id.condition]], aes(y=AcetoacetylCoA_HMGCoA, x=Cellline)) + geom_boxplot() + ggtitle("Acetoacetyl-CoA -> HMG-CoA")
ggplot(flux.list[[id.condition]], aes(y=HMGCoA_Mev, x=Cellline)) + geom_boxplot() + ggtitle("HMG-CoA -> Mevalonate")
ggplot(flux.list[[id.condition]], aes(y=HMGCRInh, x=Cellline)) + geom_boxplot() + ggtitle("HMGCR Inhibition")
ggplot(flux.list[[id.condition]], aes(y=Mev_Chol, x=Cellline)) + geom_boxplot() + ggtitle("Mevalonate -> CholER")
ggplot(flux.list[[id.condition]], aes(y=Choluptake, x=Cellline)) + geom_boxplot() + ggtitle("CholMedia -> CholER")
ggplot(flux.list[[id.condition]], aes(y=Cholshunt, x=Cellline)) + geom_boxplot() + ggtitle("CholER -> shunt")
ggplot(flux.list[[id.condition]], aes(y=AcetylCoAshunt, x=Cellline)) + geom_boxplot() + ggtitle("Acetyl-CoA -> shunt")
dev.off()
