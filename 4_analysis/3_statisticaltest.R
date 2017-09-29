library(reshape2)
library(lsr)

# select 10% best data
good.ids <- c()
for(i in unique(df.XBEST$TYPE)){
  temp=FBEST[df.XBEST$TYPE==i];
  good=which(df.XBEST$TYPE==i & FBEST <= quantile(temp, 0.5, names=FALSE)) 
  good.ids <- c(good.ids, good)
}

df.XBEST.m <- melt(df.XBEST[good.ids,], id=c("TYPE"))
df.XBEST.m <- melt(df.XBEST, id=c("TYPE"))
head(df.XBEST.m)

t.test.table <- data.frame(variable = NA, comparison1 = NA, comparison2 = NA, mean1 = NA, mean2 = NA, p.val = NA, k.val = NA, effect.size = NA)
opt.parameters <- ode_parameters$parNames[ode_parameters$index_opt_pars]

# calculate significance with t-test and kruskal-wallis test and effect size using cohensD
for(i in opt.parameters){
  data.subset <- subset(df.XBEST.m, variable == i)
  
  mean.Hek <- mean(subset(data.subset, TYPE == "Hek")$value)
  mean.Hela <- mean(subset(data.subset, TYPE == "Hela")$value)
  mean.Huh7 <- mean(subset(data.subset, TYPE == "Huh7")$value)
  mean.HepG2 <- mean(subset(data.subset, TYPE == "HepG2")$value)

  p.val.t.1 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "Hela")),  paired=FALSE, var.equal = FALSE)$p.value
  p.val.t.2 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "Huh7")),  paired=FALSE, var.equal = FALSE)$p.value
  p.val.t.3 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "HepG2")),  paired=FALSE, var.equal = FALSE)$p.value
  p.val.t.4 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hela", "Huh7")),  paired=FALSE, var.equal = FALSE)$p.value
  p.val.t.5 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hela", "HepG2")),  paired=FALSE, var.equal = FALSE)$p.value
  p.val.t.6 <- t.test(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Huh7", "HepG2")),  paired=FALSE, var.equal = FALSE)$p.value
  
  k.val.t.1 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Hek", "Hela")))$p.value
  k.val.t.2 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Hek", "Huh7")))$p.value
  k.val.t.3 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Hek", "HepG2")))$p.value
  k.val.t.4 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Hela", "Huh7")))$p.value
  k.val.t.5 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Hela", "HepG2")))$p.value
  k.val.t.6 <- kruskal.test(value ~ as.factor(TYPE), data = subset(data.subset, TYPE %in% c("Huh7", "HepG2")))$p.value
  
  effect.size.1 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "Hela")))
  effect.size.2 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "Huh7")))
  effect.size.3 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hek", "HepG2")))
  effect.size.4 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hela", "Huh7")))
  effect.size.5 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Hela", "HepG2")))
  effect.size.6 <- cohensD(value ~ TYPE, data = subset(data.subset, TYPE %in% c("Huh7", "HepG2")))
  
  t.test.table <- rbind(t.test.table, data.frame(variable = rep(i, 6), 
                                 comparison1 = c("Hek", "Hek", "Hek", "Hela", "Hela", "Huh7"), 
                                 comparison2 = c("Hela", "Huh7", "HepG2", "Huh7", "HepG2", "HepG2"),
                                 mean1 = c(mean.Hek, mean.Hek, mean.Hek, mean.Hela, mean.Hela, mean.Huh7),
                                 mean2 = c(mean.Hela, mean.Huh7, mean.HepG2, mean.Huh7, mean.HepG2, mean.HepG2),
                                 p.val = c(p.val.t.1, p.val.t.2, p.val.t.3, p.val.t.4, p.val.t.5, p.val.t.6),
                                 k.val = c(k.val.t.1, k.val.t.2, k.val.t.3, k.val.t.4, k.val.t.5, k.val.t.6),
                                 effect.size = c(effect.size.1, effect.size.2, effect.size.3, effect.size.4, effect.size.5, effect.size.6)))
  
  
}
t.test.table <- t.test.table[-1,]

# adjust for multiple testing
t.test.table$p.adj <- p.adjust(t.test.table$p.val, method="BH")
t.test.table$k.adj <- p.adjust(t.test.table$k.val, method="BH")
write.table(t.test.table, paste(Sys.Date(), "Parameter_statistics.txt", sep="_"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")

#significant parameters
length(unique(subset(t.test.table, k.adj <= 0.001)[,"variable"]))
length(unique(subset(t.test.table, k.adj <= 0.001 & effect.size >= 2)[,"variable"]))
unique(subset(t.test.table, k.adj <= 0.001 & effect.size >= 2)[,"variable"])
length(unique(subset(t.test.table, k.adj <= 0.001 & comparison1 == "Hela" & comparison2 == "Huh7")[,"variable"]))

# total parameters
length(unique(t.test.table[,"variable"]))
