
options(scipen=100)
Neu <- read.table("Sample_Neu.553x2_2607_40000_iced.AB.bed.gz", header=FALSE, sep="\t", comment.char="")
colnames(Neu) <- c("chrom", "start", "stop", "compartment", "a", "b", "c", "d", "e")
Neu_size <- Neu$stop - Neu$start

NPC <- read.table("Sample_NPC.553x2_2607_40000_iced.AB.bed.gz", header=FALSE, sep="\t", comment.char="")
colnames(NPC) <- c("chrom", "start", "stop", "compartment", "a", "b", "c", "d", "e")
NPC_size <- NPC$stop - NPC$start

Glia <- read.table("Sample_Astro.553x2_2607_40000_iced.AB.bed.gz", header=FALSE, sep="\t", comment.char="")
colnames(Glia) <- c("chrom", "start", "stop", "compartment", "a", "b", "c", "d", "e")
Glia_size <- Glia$stop - Glia$start 

png("compartment_size.png", height=2500, width=1300, res=300)
boxplot(Glia_size/1000000, NPC_size/1000000, Neu_size/1000000, ylab="Compartment size (Mb)", 
        col = c("#7570b3", "#d95f02", "#1b9e77"), names = c("Astro", "NPC", "Neu"),
        outline=FALSE, ylim=c(0,30)
        )
dev.off()

# A and B
Neu_A <- Neu[Neu$compartment == 'A',]
Neu_B <- Neu[Neu$compartment == 'B',]
NPC_A <- NPC[NPC$compartment == 'A',]
NPC_B <- NPC[NPC$compartment == 'B',]
Glia_A <- Glia[Glia$compartment == 'A',]
Glia_B <- Glia[Glia$compartment == 'B',]

Neu_A_size <- Neu_A$stop - Neu_A$start
Neu_B_size <- Neu_B$stop - Neu_B$start
NPC_A_size <- NPC_A$stop - NPC_A$start
NPC_B_size <- NPC_B$stop - NPC_B$start
Glia_A_size <- Glia_A$stop - Glia_A$start
Glia_B_size <- Glia_B$stop - Glia_B$start

lo <- read.table("overlap_loops_compartments.txt", header=TRUE, sep="\t", row.names=1)
# Remove total
rlo <- lo[,1:4]
tlo <- as.table(t(rlo))

png("compartment_loops.png", height=2500, width=1800, res=300)
barplot(tlo, col=c("red", "blue", "violet", "gray"), legend=rownames(tlo), 
        ylab= "Number of loops", main = "Compartments overlapping loop anchors")
dev.off()

# Correct for size of compartments
Neu_A_correct <- (lo["Neu", "AA"]*1000000)/(sum(Neu_A_size))
Neu_B_correct <- (lo["Neu", "BB"]*1000000)/(sum(Neu_B_size))
NPC_A_correct <- (lo["NPC", "AA"]*1000000)/(sum(NPC_A_size))
NPC_B_correct <- (lo["NPC", "BB"]*1000000)/(sum(NPC_B_size))
Glia_A_correct <- (lo["Astro", "AA"]*1000000)/(sum(Glia_A_size))
Glia_B_correct <- (lo["Astro", "BB"]*1000000)/(sum(Glia_B_size))

Neu_correct <- c(Neu_A_correct, Neu_B_correct)
NPC_correct <- c(NPC_A_correct, NPC_B_correct)
Glia_correct <- c(Glia_A_correct, Glia_B_correct)

correct_mat <- as.matrix(data.frame(Glia_correct, NPC_correct, Neu_correct, row.names=c("A", "B")))
colnames(correct_mat) <- c("Astro", "NPC", "Neu")

png("compartment_loops_size_correct.png", height=2500, width=1800, res=300)
barplot(correct_mat, col=c("red", "blue"), legend=rownames(correct_mat), 
        ylab= "Number of loops per 1 Mb for compartment", ylim=c(0,10),
        main = "Compartments overlapping loop anchors \n(compartment size corrected)")
dev.off()


