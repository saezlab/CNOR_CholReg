###################################################
# Calculate increase of cholesterol uptake capacity
###################################################
incr <- data.frame(Cellline = flux.table[flux.table$Condition == 4,"Cellline"],
                   Increase10 = flux.table[flux.table$Condition == 4,"Choluptake"] / flux.table[flux.table$Condition == 9,"Choluptake"],
                   Increase2 = flux.table[flux.table$Condition == 8,"Choluptake"] / flux.table[flux.table$Condition == 9,"Choluptake"],
                   Increase5 = flux.table[flux.table$Condition == 11,"Choluptake"] / flux.table[flux.table$Condition == 10,"Choluptake"],
                   Increase1 = flux.table[flux.table$Condition == 12,"Choluptake"] / flux.table[flux.table$Condition == 10,"Choluptake"])
incr.m <- melt(incr, id.vars = "Cellline")
incr.m$variable <- factor(incr.m$variable, levels=c("Increase2", "Increase10", "Increase1", "Increase5"))

pdf(paste(Sys.Date(), "Cholesterol_uptake_incr.pdf"))
ggplot(subset(flux.table, Condition %in% c(9) & Cellline %in% c("HEK", "Hela", "Huh7", "HepG2")), aes(y=Choluptake/(Mev_Chol+ Choluptake), x=interaction(Condition, Cellline), fill=Condition)) + geom_boxplot() + ggtitle("CholMedia -> CholER")
ggplot(subset(incr.m, Cellline %in% c("HEK", "Hela", "Huh7", "HepG2")), aes(y=value, x=interaction(variable, Cellline), fill =Cellline)) + geom_boxplot() + ggtitle("CholMedia -> CholER") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

t.test(flux.table[flux.table$Cellline == "HEK", "Choluptake"]/(flux.table[flux.table$Cellline == "HEK", "Mev_Chol"] + flux.table[flux.table$Cellline == "HEK", "Choluptake"]),
       flux.table[flux.table$Cellline == "Hela", "Choluptake"]/(flux.table[flux.table$Cellline == "Hela", "Mev_Chol"] + flux.table[flux.table$Cellline == "Hela", "Choluptake"]))
t.test(flux.table[flux.table$Cellline == "HEK", "Choluptake"]/(flux.table[flux.table$Cellline == "HEK", "Mev_Chol"] + flux.table[flux.table$Cellline == "HEK", "Choluptake"]),
       flux.table[flux.table$Cellline == "Huh7", "Choluptake"]/(flux.table[flux.table$Cellline == "Huh7", "Mev_Chol"] + flux.table[flux.table$Cellline == "Huh7", "Choluptake"]))
t.test(flux.table[flux.table$Cellline == "HEK", "Choluptake"]/(flux.table[flux.table$Cellline == "HEK", "Mev_Chol"] + flux.table[flux.table$Cellline == "HEK", "Choluptake"]),
       flux.table[flux.table$Cellline == "HepG2", "Choluptake"]/(flux.table[flux.table$Cellline == "HepG2", "Mev_Chol"] + flux.table[flux.table$Cellline == "HepG2", "Choluptake"]))
t.test(flux.table[flux.table$Cellline == "Hela", "Choluptake"]/(flux.table[flux.table$Cellline == "Hela", "Mev_Chol"] + flux.table[flux.table$Cellline == "Hela", "Choluptake"]),
       flux.table[flux.table$Cellline == "Huh7", "Choluptake"]/(flux.table[flux.table$Cellline == "Huh7", "Mev_Chol"] + flux.table[flux.table$Cellline == "Huh7", "Choluptake"]))
t.test(flux.table[flux.table$Cellline == "Hela", "Choluptake"]/(flux.table[flux.table$Cellline == "Hela", "Mev_Chol"] + flux.table[flux.table$Cellline == "Hela", "Choluptake"]),
       flux.table[flux.table$Cellline == "HepG2", "Choluptake"]/(flux.table[flux.table$Cellline == "HepG2", "Mev_Chol"] + flux.table[flux.table$Cellline == "HepG2", "Choluptake"]))
t.test(flux.table[flux.table$Cellline == "Huh7", "Choluptake"]/(flux.table[flux.table$Cellline == "Huh7", "Mev_Chol"] + flux.table[flux.table$Cellline == "Huh7", "Choluptake"]),
       flux.table[flux.table$Cellline == "HepG2", "Choluptake"]/(flux.table[flux.table$Cellline == "HepG2", "Mev_Chol"] + flux.table[flux.table$Cellline == "HepG2", "Choluptake"]))

