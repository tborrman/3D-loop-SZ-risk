df <- read.table("overlap_loops_eQTLs.txt", header=TRUE, sep="\t", comment.char="")

all_eQTL <- sum(df$GENE!= "")
all <- c(all_eQTL, length(df$GENE) - all_eQTL)

GM_eQTL <- sum(df$GMfdrDonut < -1 & df$GENE != "")
GM <- c(GM_eQTL, sum(df$GMfdrDonut < -1) - GM_eQTL)

Astro_eQTL <- sum(df$AstrofdrDonut < -1 & df$GENE != "")
Astro <- c(Astro_eQTL, sum(df$AstrofdrDonut < -1) - Astro_eQTL)

NPC_eQTL <- sum(df$NPCfdrDonut < -1 & df$GENE != "")
NPC <- c(NPC_eQTL, sum(df$NPCfdrDonut < -1) - NPC_eQTL)

Neu_eQTL <- sum(df$NeufdrDonut < -1 & df$GENE != "")
Neu <- c(Neu_eQTL, sum(df$NeufdrDonut < -1) - Neu_eQTL)

x = data.frame(all, GM, Astro, NPC, Neu, row.names= c("eQTL_loops", "non-eQTL_loops"))
x = as.matrix(x)
png("overlap_loops_eQTLs.png", height=2200, width=1700, res=300)
barplot(x, col=c("mediumorchid4", "orange"), legend=rownames(x))
text(0.75, 400, paste(round(x["eQTL_loops", "all"]/sum(x[,"all"]), 2), '%', sep=""), col="white")
text(0.75 + 1.2, 400, paste(round(x["eQTL_loops", "GM"]/sum(x[,"GM"]), 2), '%', sep=""), col="white")
text(0.75 + 1.2*2, 400, paste(round(x["eQTL_loops", "Astro"]/sum(x[,"Astro"]), 2), '%', sep=""), col="white")
text(0.75 + 1.2*3, 400, paste(round(x["eQTL_loops", "NPC"]/sum(x[,"NPC"]), 2), '%', sep=""), col="white")
text(0.75 + 1.2*4, 400, paste(round(x["eQTL_loops", "Neu"]/sum(x[,"Neu"]), 2), '%', sep=""), col="white")
dev.off()
