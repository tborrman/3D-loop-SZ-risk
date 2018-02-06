source("co-regulation_significance.R")
############################
####permutation#############
# cos <- read.table("COSgenesPositionsMatrix.txt", sep="\t", header=FALSE, comment.char="")
# cos <- cos[5:ncol(cos)]

############################################################
# Read data tables
NPC_intxns <- read.table("S12.NPC-specific_PGC_intrxns.txt", sep="\t", header=TRUE, comment.char="")
chromSizes <- read.table("chrom.sizes", sep="\t", header=FALSE)
RPKM_df <- read.table("COSgenesPositionsMatrix.txt", sep="\t", header=FALSE, comment.char="")
colnames(RPKM_df)[1:4] <- c("chrom", "start", "end", "gene")
# Separate RPKM table by chrom
RPKM_list <- split(RPKM_df, f=RPKM_df$chrom)

# Remove M chrom
chromSizes <- chromSizes[1:nrow(chromSizes) -1, ]
colnames(chromSizes) <- c("chrom", "bp")
# generate cumulative genome table
chrom <- as.character(chromSizes$chrom)
bp <- cumsum(as.numeric(chromSizes$bp))
cmChromSizes <- data.frame(chrom, bp)
###########################################################

# Distance equivalent null distribution
#NPC
count = 0
rndmeans = c()
for ( i in 1:10)
{
  #rndcos=cos[sample(nrow(cos),308),]
  rndcos <- sample_null_distribution(NPC_intxns, chromSizes, cmChromSizes, RPKM_list, RPKM_df, 308)
  rndcos <- rndcos[5:ncol(rndcos)]
  rndcos=t(rndcos)
  # Remove genes with equivalent RPKM across all samples
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  rndmeans <- c(rndmeans, rndmean)
  #print(rndmean)
  #print(n)
  if (rndmean > 0.245476) count = count+1 #0.245476 = correlation result of NPC list
}

#p = count/1000000
p = count/10

print("NPC")
print(count)
print(p)


stop()

# Random genes


#NPC
count = 0

#for ( i in 1:1000000 )
for ( i in 1:100 )
{
  rndcos=cos[sample(nrow(cos),308),]
  rndcos=t(rndcos)
  # Remove genes with equivalent RPKM across all samples
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  #print(rndmean)
  #print(n)
  if (rndmean > 0.245476) count = count+1 #0.245476 = correlation result of NPC list
}

#p = count/1000000
p = count/10000

print("NPC")
print(count)
print(p)

# Distance equivalent null distribution
x <- sample_null_distribution(NPC_intxns, chromSizes, cmChromSizes, RPKM_df, 3)


stop()



















#NEU
count = 0

for ( i in 1:1000000 )
{
  rndcos=cos[sample(nrow(cos),358),]
  
  rndcos=t(rndcos)
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  #print(rndmean)
  #print(n)
  if (rndmean > 0.2329368) count = count+1
}

p = count/1000000

print("NEU")
print(count)
print(p)




#NPC string
count = 0

for ( i in 1:1000000 )
{
  rndcos=cos[sample(nrow(cos),77),]
  
  rndcos=t(rndcos)
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  print(rndmean)
  print(n)
  if (rndmean > 0.2962837) count = count+1
}

p = count/1000000

print("NPC string")
print(count)
print(p)


#NEU string
count = 0

for ( i in 1:1000000 )
{
  rndcos=cos[sample(nrow(cos),73),]
  
  rndcos=t(rndcos)
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  #print(rndmean)
  if (rndmean > 0.2877246) count = count+1
}

p = count/1000000

print("NEU string")
print(count)
print(p)



#NPC string medium
count = 0

for ( i in 1:1000000 )
{
  rndcos=cos[sample(nrow(cos),139),]
  
  rndcos=t(rndcos)
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  print(rndmean)
  print(n)
  if (rndmean > 0.2725747) count = count+1
}

p = count/1000000

print("NPC string medium")
print(count)
print(p)



#NEU string medium
count = 0

for ( i in 1:1000000 )
{
  rndcos=cos[sample(nrow(cos),138),]
  
  rndcos=t(rndcos)
  zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
  rndcos <- rndcos[, !zv]
  n=length(colnames(rndcos))
  
  crnd <- cor(rndcos[,1:n],use="complete.obs")
  #pheatmap(crnd, color = colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(100), border_color=NA, breaks=myBreaks)
  
  crnd1=abs(crnd)
  rndmean=mean(crnd1)
  #print(rndmean)
  #print(n)
  if (rndmean > 0.2654931) count = count+1
  #print(count)
}

p = count/1000000

print("NEU string medium")
print(count)
print(p)






####################
#######graph########

geneno = 77 #number of genes in tested gene list
histobase=replicate(10000,
                    {
                      rndcos=cos[sample(nrow(cos),geneno),] 
                      rndcos=t(rndcos)
                      zv <- apply(rndcos, 2, function(x) length(unique(x)) == 1)
                      rndcos <- rndcos[, !zv]
                      n=length(colnames(rndcos))
                      crnd <- cor(rndcos[,1:n],use="complete.obs")
                      crnd=abs(crnd)
                      return(mean(crnd))
                    })

hist(histobase,breaks=20,freq=FALSE,xlim=c(0.12,0.31), xlab="absolute correlation coefficient", col ="lightgrey")
abline(v=0.2963, col="red", lwd=2)
