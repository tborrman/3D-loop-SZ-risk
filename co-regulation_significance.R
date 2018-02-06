
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



is_cis <- function(row) {
  # Determine if significant interaction is in
  # cis or trans
  # Returns: boolean; True if in cis"
  chr1 <- sub(":.*", "", row$anchor.bin.coord)
  chr2 <- sub(":.*", "", row$target.bin.coord)
  if(chr1 == chr2){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

get_distance <- function(row) {
  # Returns: bp distance between significant interaction
  start <- as.numeric(sub(".*:(.*)-.*", "\\1", row$anchor.bin.coord))
  end <- as.numeric(sub(".*:(.*)-.*", "\\1", row$target.bin.coord))
  d <- abs(start - end)
  return(d)
}

rand_genomic_pos <- function(d, c, cm) {
  # Args: '
  # d = bp distance
  # c = chromosome size table
  # cm = cumulative chromosome size table 
  # Returns: list coords
  #   chrom = chromosome
  #   g1 = random genomic coordinate
  #   g2 = paired genomic coordinate d bp from g1
  total <- cm$bp[length(cm$bp)]
  chrStarts <- c(1,cm$bp[1:length(cm$bp) -1] + 1)
  rand_pos <- sample(total, 1)
  chr_idx <- max(which(chrStarts <= rand_pos))
  chrom <- as.character(c$chrom[chr_idx])
  g1 <- rand_pos - chrStarts[chr_idx] + 1
  g2 <- g1 + d
  coords <- list(chrom, g1, g2)
  if (g2 > c$bp[chr_idx]) {
    print(coords)
    coords<- rand_genomic_pos(d, c, cm)
  }
  return(coords)
}

overlap <- function(a, b){
  # Return boolean for whether interval b
  # overlaps interval a
  
  if (((a[1] <= b[1]) & (b[1]<= a[2])) | ((b[1] <= a[1]) & (a[1]<= b[2]))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_overlap_gene_idxs <- function(RPKM_df, rand_coord) {
  # Return row numbers of RPKM_df for genes
  # overlapping input random coordinates
  # Args:
  # RPKM_df = RPKM table for genes
  # rand_coord = list of [chrom, coord1, coord2]
  # Returns: 
  # vector idx = row indexes for genes that overlap 
  # random genome coordinates
  rand_chrom <- rand_coord[[1]]
  # 10kb bins
  rand_g1 <- c(rand_coord[[2]], rand_coord[[2]] + 10000)
  rand_g2 <- c(rand_coord[[3]], rand_coord[[3]] + 10000)
  #print(rand_coord)
  #print(rand_g1)
  #print(rand_g2)
  idx = c()
  chrom_gene_df <-  RPKM_list[[rand_chrom]]
  for (i in 1:nrow(chrom_gene_df)) {
    gene_pos <- c(chrom_gene_df[i,]$start, chrom_gene_df[i,]$end)
    # if (i %% 100 == 0) {
    #   print(i)
    # }
    if (overlap(gene_pos, rand_g1) | overlap(gene_pos, rand_g2)) {
      # print(rand_chrom)
      # print(rand_g1)
      # print(rand_g2)
      # print(chrom_gene_df[i,]$chrom)
      # print(gene_pos)
      # print(chrom_gene_df[i,1:5])
      # print(rownames(chrom_gene_df[i,]))
      idx <- c(idx, as.numeric(rownames(chrom_gene_df[i,])))
    }
  }
  return(idx)
}



# Get random PGC interacion
r <- sample(1:nrow(NPC_intxns), 1)
randRow = NPC_intxns[r,]
if (is_cis(randRow)) {
  d <- get_distance(randRow)
  rand_coord <- rand_genomic_pos(d, chromSizes, cmChromSizes)
  print(rand_coord)
  rand_idxs <- get_overlap_gene_idxs(RPKM_df, rand_coord)
  print(rand_idxs)

} else {
  print("ERROR: trans interaction")
  stop()
}




