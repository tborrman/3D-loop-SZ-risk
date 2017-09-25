df <- read.table("gene_table", sep="\t", header=TRUE, quote="", comment.char="")
NPC_df <- read.table("NPC_overlapping_genes.txt", sep="\t", header=FALSE, quote="", comment.char="")
tss_df <- read.table("NPC_overlapping_tss.txt", sep="\t", header=FALSE, quote="", comment.char="")
full_df <- read.table("NPC_overlapping_genes_full_loop.txt", sep="\t", header=FALSE, quote="", comment.char="")
colnames(tss_df) <- colnames(df)
colnames(NPC_df) <- colnames(df)
colnames(full_df) <- colnames(df)

z.test = function(x,mu,popvar){
  z.score <- (mean(x)-mu)/(popvar/sqrt(length(x)))
  print(z.score)
  one.tail.p <- pnorm(abs(z.score),lower.tail = FALSE)
  return(one.tail.p)
}

pop_expr = log10(df$X2607.1.AN.2 + 1)
loop_expr = log10(NPC_df$X2607.1.AN.2 + 1)

p <- z.test(loop_expr, mean(pop_expr), sd(pop_expr))
print(p)


tss_expr = log10(tss_df$X2607.1.AN.2 + 1)

p2 <-z.test(tss_expr, mean(pop_expr), sd(pop_expr))
print(p2)

full_expr = log10(full_df$X2607.1.AN.2 + 1)

p3 <-z.test(full_expr, mean(pop_expr), sd(pop_expr))
print(p3)

x <- sample(1:length(pop_expr), length(loop_expr))
y <- sample(1:length(pop_expr), length(loop_expr))
z <- sample(1:length(pop_expr), length(loop_expr))

p4 <-z.test(pop_expr[x], mean(pop_expr), sd(pop_expr))
print(p4)

png('FPKM_boxplots.png', height=2500, width=1300, res=300)
boxplot(pop_expr, full_expr, loop_expr, tss_expr, names=c("All", "Loop", "Anchor", "TSS"), 
        ylab=expression('log'[10]*'(FPKM + 1)'), xlab='Gene sets',
        col=c("gray","darkslategrey", "#1b9e77", "darkseagreen3"))
dev.off()
